using NetworkDynamics
using Aqua
using ExplicitImports
using Test

@testset "Package Quality Tests" begin
    # print_explicit_imports(NetworkDynamics)
    @test check_no_implicit_imports(NetworkDynamics) === nothing
    @test check_no_stale_explicit_imports(NetworkDynamics, ignore=(:Symbolics,)) === nothing
    Aqua.test_all(NetworkDynamics;
        ambiguities=false,
        # Hungarian and Moshi are only needed for the MTK extension (Moshi used to be
        # pulled in transitively via SciMLBase<3, which now makes it an optional ext).
        stale_deps=(; ignore=[:Hungarian, :Moshi]),
        persistent_tasks=false)
    @test_broken isempty(Docs.undocumented_names(NetworkDynamics))
end
