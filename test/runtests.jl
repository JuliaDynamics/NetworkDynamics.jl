using NetworkDynamics
using Test
using Graphs
using SafeTestsets

using NetworkDynamics: VertexBatch, parameter_range

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")


@testset "ND Prototype Tests" begin

@testset "Test Component Library" begin
    using NetworkDynamics: compf
    a = Lib.diffusion_edge()
    b = Lib.diffusion_edge()
    @test compf(a) == compf(b)
    a = Lib.diffusion_edge_closure()
    b = Lib.diffusion_edge_closure()
    @test compf(a) != compf(b)

    fid = Lib.diffusion_edge_fid()
    ode = Lib.diffusion_odeedge()

    odevert = Lib.diffusion_vertex()
    odevert_const = Lib.diffusion_vertex_constraint()

    kura_edge   = Lib.kuramoto_edge()
    kura_second = Lib.kuramoto_second()
end

include("utils_test.jl")
include("construction_test.jl")

@safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
@safetestset "Symbolic Indexing Tests" begin include("symbolicindexing_test.jl") end

end
