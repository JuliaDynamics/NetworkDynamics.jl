using NetworkDynamics, Graphs
using Chairmarks: @b

using NetworkDynamics: StateBufIdx, OutBufIdx, ExtMap, collect_externals!
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as Dt
using SciCompDSL
using Random
using CUDA, Adapt

@testset "test extmap performance" begin
    N = 10_000
    u = rand(N)
    o = rand(N)
    _map = Union{StateBufIdx,OutBufIdx}[rand((StateBufIdx(i), OutBufIdx(i))) for i in 1:N]
    map = ExtMap(Random.shuffle(_map))

    exbuf = zeros(N)
    fill!(exbuf, NaN)
    b = @b collect_externals!($map, $exbuf, $u, $o) # 6.583μs
    @test iszero(b.allocs)

    if CUDA.functional()
        for to in (CuArray{Float64}, CuArray{Float32})
            map_d = adapt(to, map)
            exbuf_d = adapt(to, zeros(N))
            u_d = adapt(to, u)
            o_d = adapt(to, o)
            b = @b collect_externals!($map_d, $exbuf_d, $u_d, $o_d)
            display(b)
            @test Vector(exbuf_d) ≈ exbuf
        end
    end
end

@testset "test construction of elements with external inputs" begin
    function fvext(dv, v, ein, ext, p, t)
        dv[1] = ext[1]
    end
    function gvext(out, v, ein, ext, p, t)
        out .= ext
    end
    v_noff = VertexModel(f=fvext, g=1, dim=1, extin=[VIndex(2,:a)])
    v_ff   = VertexModel(f=fvext, g=gvext, outdim=1, dim=1, extin=[VIndex(2,:a)])

    function feext(de, e, vsrc, vdst, ext, p ,t)
        de[1] = ext[1]
    end
    function geext(outsrc, outdst, e, vsrc, vdst, ext, p, t)
        outsrc .= ext
        outdst .= ext
    end
    e_noff = EdgeModel(f=feext, g=AntiSymmetric(1), dim=1, extin=[EIndex(2,:a)])
    e_ff   = EdgeModel(f=feext, g=geext, outdim=1, dim=1, extin=[EIndex(2,:a)])

    @mtkmodel EdgeExt begin
        @variables begin
            srcin(t)
            dstin(t)
            extin1(t)
            extin2(t)
            out(t)
        end
        @equations begin
            out ~ extin1 + extin2
        end
    end
    @named edge = EdgeExt()
    edgem = EdgeModel(edge, :srcin, :dstin, Directed(:out); extin=[:extin1 => VIndex(2,:a), :extin2 => VIndex(2,:a)])
end

@testset "external input referencing an output which is not a state" begin
    # vertex 2 has state x and output o = 2x, which is not a state
    vsrc = VertexModel(f=(dv, v, esum, p, t) -> dv[1] = 0, g=(out, v, p, t) -> out[1] = 2v[1],
                       sym=[:x], outsym=[:o], dim=1, outdim=1)
    # vertex 1 integrates o, the MTK edge outputs o
    vdst = VertexModel(f=(dv, v, esum, ext, p, t) -> dv[1] = ext[1], g=1, sym=[:y], dim=1,
                       extin=[VIndex(2, :o)])
    @mtkmodel EdgeOutExt begin
        @variables begin
            srcin(t)
            dstin(t)
            ext(t)
            out(t)
        end
        @equations begin
            out ~ ext
        end
    end
    @named edge = EdgeOutExt()
    em = EdgeModel(edge, :srcin, :dstin, AntiSymmetric(:out); extin=[:ext => VIndex(2, :o)])

    nw = Network(path_graph(2), [vdst, vsrc], em)
    s = NWState(nw)
    s.v[1, :y] = 0.0
    s.v[2, :x] = 1.5
    du = zeros(dim(nw))
    nw(du, uflat(s), pflat(s), 0.0)
    @test s.v[2, :o] == 3.0
    @test NWState(nw, du).v[1, :y] == 3.0
    @test s.e[1, :out] == 3.0
end
