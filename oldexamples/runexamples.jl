begin
    using Pkg
    Pkg.activate(@__DIR__)
    println("----------------------Testing ND Examples----------------------")
    for file in
        ("accessing_edge_values.jl",
         "cascading_failure.jl",
         "diffusion.jl",
         "directed_and_weighted_graphs.jl",
         "examples_without_networkdynamics.jl",
         "fixpoint.jl",
         "getting_started_with_DDEs.jl",
         "getting_started.jl",
         "kuramoto_inertial.jl",
         "kuramoto_oscillators.jl",
         "kuramoto_plasticity.jl",
         "power_network_example.jl",
         "sde.jl")


        println("\nTesting ", file, ".\n")
        begin
            try
                include(file)
            catch e
                println(e)
                print("!!! ERROR in ", file, ". Continuing with the next example...\n")
            end
        end
    end
end
