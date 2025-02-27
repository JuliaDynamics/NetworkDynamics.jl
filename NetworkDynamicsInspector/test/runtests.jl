using Test
using NetworkDynamics
using NetworkDynamicsInspector
using ExplicitImports
using Aqua

@testset "Package Quality Tests" begin
    @test check_no_implicit_imports(NetworkDynamicsInspector) === nothing
    @test check_no_stale_explicit_imports(NetworkDynamicsInspector) === nothing

    Aqua.test_all(NetworkDynamicsInspector;
        ambiguities=false,
        stale_deps=true,
        persistent_tasks=false,
        piracies = (;treat_as_own = [NetworkDynamics.extract_nw])
    )

    @test_broken isempty(Docs.undocumented_names(NetworkDynamicsInspector))
end
