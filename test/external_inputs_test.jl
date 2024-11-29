using NetworkDynamics, Graphs
using Chairmarks: @b

using NetworkDynamics: StateBufIdx, OutBufIdx, ExtMap, collect_externals!
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
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
