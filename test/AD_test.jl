using NetworkDynamics
using Graphs
using Test
using DifferentiationInterface
using DifferentiationInterfaceTest
using SparseConnectivityTracer

import ForwardDiff, FiniteDiff, ReverseDiff, Enzyme
# import Mooncake
import Enzyme: EnzymeCore

@__MODULE__()==Main ? includet(joinpath(pkgdir(NetworkDynamics), "test", "ComponentLibrary.jl")) : (const Lib = Main.Lib)

g = complete_graph(4)
vf = [Lib.kuramoto_second(), Lib.diffusion_vertex(), Lib.kuramoto_second(), Lib.diffusion_vertex()]
ef = [Lib.diffusion_odeedge(),
      Lib.kuramoto_edge(),
      Lib.kuramoto_edge(),
      Lib.diffusion_edge_fid(),
      Lib.diffusion_odeedge(),
      Lib.diffusion_edge_fid()]
nw = Network(g, vf, ef)

x0 = rand(dim(nw))
p0 = NWParameter(nw)
p0.e[2:3,:K] .= 1.0
pf = pflat(p0)
@assert !any(isnan, pf) "test setup broken: NaN parameters give NaN jacobians"

# derivative w.r.t. state and w.r.t. parameters, in both the allocating and the
# mutating form
fx     = x -> (dx = similar(x); nw(dx, x, pf, 0.0); dx)
fp     = p -> (dx = similar(p, length(x0)); nw(dx, x0, p, 0.0); dx)
fx!    = (dx, x) -> nw(dx, x, pf, 0.0)
fp!    = (dx, p) -> nw(dx, x0, p, 0.0)

Jx_ref = jacobian(fx, AutoFiniteDiff(), x0)
Jp_ref = jacobian(fp, AutoFiniteDiff(), pf)

@testset "AD backends on the network rhs" begin
    # Enzyme is listed without `set_runtime_activity`: the rhs must stay free of
    # mixed-activity captures, and a regression there shows up as an
    # EnzymeRuntimeActivityError here.
    backends = [
        "ForwardDiff" => AutoForwardDiff(),
        "ReverseDiff" => AutoReverseDiff(),
        # "Mooncake"    => AutoMooncake(),
        "Enzyme fwd"  => AutoEnzyme(; mode=EnzymeCore.Forward, function_annotation=EnzymeCore.Const),
        "Enzyme rev"  => AutoEnzyme(; mode=EnzymeCore.Reverse, function_annotation=EnzymeCore.Const),
    ]
    @testset "$name" for (name, backend) in backends
        @testset "d/dx out-of-place" begin
            @test jacobian(fx, backend, x0) ≈ Jx_ref rtol=1e-5
        end
        @testset "d/dx in-place" begin
            @test jacobian(fx!, zeros(dim(nw)), backend, x0) ≈ Jx_ref rtol=1e-5
        end
        @testset "d/dp in-place" begin
            @test jacobian(fp!, zeros(dim(nw)), backend, pf) ≈ Jp_ref rtol=1e-5
        end
        @testset "d/dp out-of-place" begin
            if name == "Enzyme rev"
                # `similar(p, ...)` inside the closure trips Enzyme reverse mode with
                # "Conversion of boxed type Vector{Float64} is not allowed"
                @test_broken try
                    jacobian(fp, backend, pf) ≈ Jp_ref
                catch
                    false
                end
            else
                @test jacobian(fp, backend, pf) ≈ Jp_ref rtol=1e-5
            end
        end
    end
end

# @testset "DifferentiationInterfaceTest scenarios" begin
#     scenarios = [Scenario{:jacobian, :in}(fx, x0; res1=Jx_ref),
#                  Scenario{:jacobian, :in}(fp, pf; res1=Jp_ref)]
#     backends = [AutoForwardDiff(), AutoReverseDiff(), AutoMooncake(),
#                 AutoEnzyme(; mode=EnzymeCore.Forward, function_annotation=EnzymeCore.Const)]
#     test_differentiation(
#         backends,
#         scenarios,
#         correctness=true,
#         type_stability=:none,
#         detailed=true,
#         benchmark=:none,
#     )
# end

@testset "sparsity tracer" begin
    detector = TracerSparsityDetector()
    @test jacobian_sparsity(fx, x0, detector) isa AbstractMatrix
    @test jacobian_sparsity(fp, pf, detector) isa AbstractMatrix
end
