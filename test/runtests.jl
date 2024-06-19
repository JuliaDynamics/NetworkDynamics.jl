using Test
using SafeTestsets
using Pkg
using NetworkDynamics
using SciMLBase
using InteractiveUtils

"""
Test utility, which rebuilds the Network with all different execution styles and compares the
results of the coreloop.
"""
function test_execution_styles(prob)
    styles = [KAExecution{true}(), KAExecution{false}(), SequentialExecution{true}(), SequentialExecution{false}()]
    # styles = [KAExecution{false}(), SequentialExecution{true}()]
    unmatchedstyles = filter(subtypes(NetworkDynamics.ExecutionStyle)) do abstractstyle
        !any(s -> s isa abstractstyle, styles)
    end
    @assert isempty(unmatchedstyles) "Some ExecutionStyle won't be tested: $unmatchedstyles"

    aggregators = [NaiveAggregator, NNlibScatter, KAAggregator, SequentialAggregator, PolyesterAggregator]
    unmatchedaggregators = filter(subtypes(NetworkDynamics.Aggregator)) do abstractaggregator
        !any(s -> s <: abstractaggregator, aggregators)
    end
    @assert isempty(unmatchedaggregators) "Some AggrgationStyle won't be tested: $unmatchedaggregators"

    @assert prob isa ODEProblem "test_execution_styles only works for ODEProblems"

    u = copy(prob.u0)
    du = zeros(eltype(u), length(u))
    t = 0.0
    p = copy(prob.p)
    nw = prob.f.f
    nw(du, u, p, t)
    @test u==prob.u0
    @test p==prob.p

    @testset "Test Execution Styles and Aggregators" begin
        for execution in styles
            for aggregator in aggregators
                _nw = Network(nw; execution, aggregator=aggregator(nw.layer.aggregator.f))
                _du = zeros(eltype(u), length(u))
                try
                    _nw(_du, u, p, t)
                catch e
                    println("Error in $execution with $aggregator: $e")
                    @test false
                    continue
                end
                issame = _du ≈ du

                if !issame
                    println("$execution with $aggregator lead to different results: extrema(Δ) = $(extrema(_du - u))")
                end
                @test issame
            end
        end
    end
end

@testset "NetworkDynamics Tests" begin
    @safetestset "utils test" begin include("utils_test.jl") end
    @safetestset "construction test" begin include("construction_test.jl") end
    @safetestset "Aggregation Tests" begin include("aggregators_test.jl") end
    @safetestset "Symbolic Indexing Tests" begin include("symbolicindexing_test.jl") end

    @safetestset "Diffusion test" begin include("diffusion_test.jl") end
    @safetestset "inhomogeneous test" begin include("inhomogeneous_test.jl") end
    @safetestset "massmatrix test" begin include("massmatrix_test.jl") end
    @safetestset "doctor test" begin include("doctor_test.jl") end
    @safetestset "initialization test" begin include("initialization_test.jl") end
end

@testset "Test Doc Examples" begin
    @info "Activate doc environment and test examples"
    Pkg.activate(joinpath(pkgdir(NetworkDynamics), "docs"))
    Pkg.develop(path=pkgdir(NetworkDynamics))
    Pkg.instantiate()

    examples = joinpath(pkgdir(NetworkDynamics), "docs", "examples")
    for file in readdir(examples; join=true)
        endswith(file, ".jl") || continue
        name = basename(file)

        @info "Test $name"
        if name == "kuramoto_delay.jl"
            @test_broken false
            continue
        end
        eval(:(@safetestset $name begin include($file) end))
    end
end
