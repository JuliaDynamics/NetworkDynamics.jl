using NetworkDynamics, Graphs
using Chairmarks: @b

using NetworkDynamics: StateBufIdx, OutBufIdx, ExtMap, collect_externals!
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

@testset "test construction of elements with external inputs" begin
    function fext(dv, v, ein, p, t, ext)
        dv[1] = ext[1]
    end
    c = VertexModel(f=fext, g=1, dim=1, extin=[VIndex(18,:a)])
end
