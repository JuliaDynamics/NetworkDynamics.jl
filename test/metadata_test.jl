using NetworkDynamics
using DataFrames
using Graphs

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

function basenetwork()
    g = SimpleGraph([0 1 1 0 1;
                     1 0 1 1 0;
                     1 1 0 1 0;
                     0 1 1 0 1;
                     1 0 0 1 0])

    vs = [Lib.swing_mtk() for _ in 1:5];
    set_default!(vs[1], :Pmech, -1)
    set_default!(vs[2], :Pmech, 1.5)
    set_default!(vs[3], :Pmech, -1)
    set_default!(vs[4], :Pmech, -1)
    set_default!(vs[5], :Pmech, 1.5)

    ls = [Lib.line_mtk() for _ in 1:7];
    nw = Network(g, vs, ls)
    sinit = NWState(nw)
    s0 = find_fixpoint(nw)
    set_defaults!(nw, s0)
    nw
end

@testset "Symbol existence checks in metadata functions" begin
    # Create a test vertex model
    vm = Lib.swing_mtk()

    # Test with existing symbols (should not throw)
    set_metadata!(vm, :θ, :test_key, "test_value")
    @test get_metadata(vm, :θ, :test_key) == "test_value"
    @test has_metadata(vm, :θ, :test_key)

    # Test with non-existing symbols (should throw ArgumentError)
    @test_throws ArgumentError set_metadata!(vm, :nonexistent, :test_key, "test_value")
    @test_throws ArgumentError get_metadata(vm, :nonexistent, :test_key)
    @test_throws ArgumentError has_metadata(vm, :nonexistent, :test_key)

    # Test with different types of symbols
    # State symbol
    set_metadata!(vm, :ω, :test_key, "test_value")
    # Parameter symbol
    set_metadata!(vm, :M, :test_key, "test_value")
    # Input symbol
    set_metadata!(vm, :P, :test_key, "test_value")
    # Output symbol
    set_metadata!(vm, :θ, :test_key, "test_value")
    # Observed symbol
    set_metadata!(vm, :Pdamping, :test_key, "test_value")

    # Test the pair version of set_metadata!
    set_metadata!(vm, :θ, :test_key => "new_value")
    @test get_metadata(vm, :θ, :test_key) == "new_value"
    @test_throws ArgumentError set_metadata!(vm, :nonexistent, :test_key => "value")

    # Test auto-generated metadata methods
    set_default!(vm, :θ, 1.0)
    @test_throws ArgumentError set_default!(vm, :nonexistent, 1.0)
    set_guess!(vm, :M, 2.0)
    @test_throws ArgumentError set_guess!(vm, :nonexistent, 2.0)
    set_init!(vm, :ω, 0.5)
    @test_throws ArgumentError set_init!(vm, :nonexistent, 0.5)
    set_bounds!(vm, :D, (0.0, 1.0))
    @test_throws ArgumentError set_bounds!(vm, :nonexistent, (0.0, 1.0))
end

@testset "Network metadata access" begin
    # Create a test network
    nw = basenetwork()

    # Test metadata access using SymbolicIndex
    vidx = VIndex(1, :θ)

    # Test setting metadata on a network element
    set_metadata!(nw, vidx, :test_key, "test_value")
    @test has_metadata(nw, vidx, :test_key)
    @test get_metadata(nw, vidx, :test_key) == "test_value"

    # Test the pair version
    set_metadata!(nw, vidx, :another_key => "another_value")
    @test get_metadata(nw, vidx, :another_key) == "another_value"

    # Test auto-generated metadata methods with network
    set_default!(nw, vidx, 2.0)
    @test has_default(nw, vidx)
    @test get_default(nw, vidx) == 2.0

    # Test combined metadata accessors
    @test NetworkDynamics.has_default_or_init(nw, vidx)
    @test NetworkDynamics.get_default_or_init(nw, vidx) == 2.0

    # Test with non-existent symbol
    invalid_idx = VIndex(1, :nonexistent)
    @test_throws ArgumentError set_metadata!(nw, invalid_idx, :test_key, "test_value")
    @test_throws ArgumentError get_metadata(nw, invalid_idx, :test_key)
    @test_throws ArgumentError has_metadata(nw, invalid_idx, :test_key)
end

@testset "Component-level metadata with network access" begin
    nw = basenetwork()

    # Test setting component metadata
    set_metadata!(nw, VIndex(1), :test_comp_key, "test_comp_value")
    @test has_metadata(nw, VIndex(1), :test_comp_key)
    @test get_metadata(nw, VIndex(1), :test_comp_key) == "test_comp_value"

    # Test graphelement functions
    set_position!(nw, VIndex(1), (0.5, 0.5))
    @test has_position(nw, VIndex(1))
    @test get_position(nw, VIndex(1)) == (0.5, 0.5)

    set_marker!(nw, VIndex(1), :circle)
    @test has_marker(nw, VIndex(1))
    @test get_marker(nw, VIndex(1)) == :circle

    # Test callbacks
    cb = ContinousComponentCallback(
        ComponentCondition(nothing, [], []),
        ComponentAffect(nothing, [], [])
    )
    set_callback!(nw, VIndex(1), cb)
    @test has_callback(nw, VIndex(1))
    @test get_callbacks(nw, VIndex(1)) == (cb,)

    cb2 = ContinousComponentCallback(
        ComponentCondition(nothing, [], []),
        ComponentAffect(nothing, [], [])
    )
    add_callback!(nw, VIndex(1), cb2)
    @test length(get_callbacks(nw, VIndex(1))) == 2

    # Test bulk metadata retrieval functions
    syms = [:θ, :ω]
    @test length(NetworkDynamics.get_defaults(nw, VIndex(1), syms)) == 2
    @test length(NetworkDynamics.get_guesses(nw, VIndex(1), syms)) == 2
    @test length(NetworkDynamics.get_defaults_or_inits(nw, VIndex(1), syms)) == 2
end

@testset "describe_ functions" begin
    nw = basenetwork()
    # Test describe_vertices
    df_vertices = describe_vertices(nw)
    @test allequal(df_vertices.batch)

    # Check DataFrame structure
    @test df_vertices isa DataFrame
    @test size(df_vertices, 1) == 5  # 5 vertices in basenetwork

    # Check essential columns are present
    @test :idx in propertynames(df_vertices)
    @test :name in propertynames(df_vertices)
    @test :batch in propertynames(df_vertices)

    # Check parameter values like Pmech appear in the DataFrame
    @test :Pmech in propertynames(df_vertices)
    @test df_vertices[1, :Pmech] == -1.0
    @test df_vertices[2, :Pmech] == 1.5
    @test df_vertices[5, :Pmech] == 1.5

    # Test describe_edges
    df_edges = describe_edges(nw)
    @test allequal(df_edges.batch)

    # Check DataFrame structure
    @test df_edges isa DataFrame
    @test size(df_edges, 1) == 7  # 7 edges in basenetwork

    # Check essential columns are present
    @test :idx in propertynames(df_edges)
    @test :srcdst in propertynames(df_edges)
    @test :name in propertynames(df_edges)
    @test :batch in propertynames(df_edges)

    # Test with parameter exclusion
    df_vertices_no_params = describe_vertices(nw; parameters=false)
    @test !(:Pmech in propertynames(df_vertices_no_params))

    # Test with state inclusion/exclusion
    df_edges_no_states = describe_edges(nw; states=false)
    for state_sym in NetworkDynamics.sym(first(nw.im.edgem))
        @test !(state_sym in propertynames(df_edges_no_states))
    end

    # Test with custom metadata
    set_metadata!(nw[VIndex(1)], :custom_key, "test_value")
    df_with_custom = describe_vertices(nw, :custom_key => vm -> has_metadata(vm, :custom_key) ? get_metadata(vm, :custom_key) : missing)
    @test :custom_key in propertynames(df_with_custom)
    @test df_with_custom[1, :custom_key] == "test_value"
    @test ismissing(df_with_custom[2, :custom_key])
end
