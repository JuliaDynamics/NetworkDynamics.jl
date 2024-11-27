using NetworkDynamics, Graphs
using Chairmarks: @b

using NetworkDynamics: StateBufIdx, OutBufIdx, ExtMap, collect_externals!
using Random

@testset "test extmap performance" begin
    N = 10_000
    u = rand(N)
    o = rand(N)
    _map = Union{StateBufIdx,OutBufIdx}[rand((StateBufIdx(i), OutBufIdx(i))) for i in 1:N]
    map = ExtMap(Random.shuffle(_map))

    exbuf = zeros(N)
    fill!(exbuf, NaN)
    b = @b collect_externals!(map, exbuf, u, o) # 6.583μs
    @test iszero(b.allocs)

    _map = [StateBufIdx(i) for i in 1:N]
    map = ExtMap(Random.shuffle(_map))
    exbuf = zeros(N)
    fill!(exbuf, NaN)
    b = @b collect_externals!(map, exbuf, u, o) # 5.6μs
    @test iszero(b.allocs)
end

@testset "test construction of elements with external inputs" begin
    function fext(dv, v, ein, p, t, ext)
        dv[1] = ext[1]
    end
    c = VertexModel(f=fext, g=1, dim=1, extsym=[VIndex(18,:a)])
end
