using NetworkDynamics
using DataFrames
using Graphs
using Test

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

@testset "Graphelement functions for edges" begin
    # Create edge model
    em = Lib.line_mtk()

    # Set graphelement on edge model directly
    set_graphelement!(em, (;src=1, dst=2))
    @test has_graphelement(em)
    @test get_graphelement(em) == (;src=1, dst=2)

    # Test using symbolically named vertices
    set_graphelement!(em, (;src=:source, dst=:target))
    @test get_graphelement(em) == (;src=:source, dst=:target)

    # Test alternative pair syntax
    set_graphelement!(em, 3 => 4)
    @test get_graphelement(em) == (;src=3, dst=4)

    # Test with Network access
    nw = basenetwork()
    eidx = EIndex(1)

    set_graphelement!(nw, eidx, 5 => 6)
    @test has_graphelement(nw, eidx)
    @test get_graphelement(nw, eidx) == (;src=5, dst=6)
end

@testset "Initial state retrieval" begin
    # Test direct ComponentModel access
    vm = Lib.swing_mtk()
    set_default!(vm, :θ, 0.5)
    set_init!(vm, :ω, 0.0)

    # Test get_initial_state for a single symbol
    @test get_initial_state(vm, :θ) == 0.5
    @test get_initial_state(vm, :ω) == 0.0

    # Test get_initial_state for multiple symbols
    states = get_initial_state(vm, [:θ, :ω])
    @test states == [0.5, 0.0]

    # Test with missing value
    @test get_initial_state(vm, :Pmech, missing_val=:missing) == :missing

    # Test with Network access
    nw = basenetwork()
    vidx = VIndex(1, :θ)
    eidx = EIndex(1, :P)
    set_default!(nw, vidx, 0.3)
    @test get_initial_state(nw, vidx) == 0.3
    set_init!(nw, eidx, 0.2)
    @test get_initial_state(nw, eidx) == 0.2
end

@testset "Aliased component detection" begin
    # Create network with aliased components
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)

    vm = Lib.swing_mtk()
    em = Lib.line_mtk()

    # Create network with same component used multiple times (aliased)
    nw = Network(g, [vm, vm, vm], [em, em])

    # Initially no changes
    @test !NetworkDynamics.aliased_changed(nw; warn=false)

    # Make a change to the aliased component
    set_default!(vm, :θ, 0.5)

    # Should detect the change
    @test NetworkDynamics.aliased_changed(nw; warn=false)

    # Test with edges
    g2 = SimpleGraph(3)
    add_edge!(g2, 1, 2)
    add_edge!(g2, 2, 3)

    em_shared = Lib.line_mtk()
    nw2 = Network(g2, [Lib.swing_mtk() for _ in 1:3], em_shared)

    # Make a change to the aliased edge component
    set_default!(em_shared, :K, 0.3)

    # Should detect the change
    @test NetworkDynamics.aliased_changed(nw2; warn=false)
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

    # Test delete_metadata! for component-wide metadata
    @test delete_metadata!(nw, VIndex(1), :test_comp_key)
    @test !has_metadata(nw, VIndex(1), :test_comp_key)
    # Test deleting non-existent metadata
    @test !delete_metadata!(nw, VIndex(1), :nonexistent_key)

    # Test graphelement functions
    set_position!(nw, VIndex(1), (0.5, 0.5))
    @test has_position(nw, VIndex(1))
    @test get_position(nw, VIndex(1)) == (0.5, 0.5)

    set_marker!(nw, VIndex(1), :circle)
    @test has_marker(nw, VIndex(1))
    @test get_marker(nw, VIndex(1)) == :circle

    # Test callbacks
    cb = ContinuousComponentCallback(
        ComponentCondition(nothing, [], []),
        ComponentAffect(nothing, [], [])
    )
    set_callback!(nw, VIndex(1), cb)
    @test has_callback(nw, VIndex(1))
    @test get_callbacks(nw, VIndex(1)) == (cb,)

    cb2 = ContinuousComponentCallback(
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

@testset "Bulk metadata retrieval functions" begin
    # Create a test ComponentModel
    vm = Lib.swing_mtk()

    # Set up some test data
    set_default!(vm, :θ, 0.5)
    set_init!(vm, :ω, 0.1)
    set_guess!(vm, :M, 5.0)
    delete_default!(vm, :ω)

    # Test get_defaults for multiple symbols
    syms = [:θ, :ω, :Pmech]
    defaults = NetworkDynamics.get_defaults(vm, syms)
    @test defaults[1] == 0.5  # θ has default
    @test defaults[2] === nothing  # ω has no default
    @test defaults[3] === nothing  # Pmech has no default

    # Test get_guesses
    syms = [:θ, :M, :D]
    guesses = NetworkDynamics.get_guesses(vm, syms)
    @test guesses[1] === nothing  # θ has no guess
    @test guesses[2] == 5.0  # M has guess
    @test guesses[3] === nothing  # D has no guess

    # Test get_defaults_or_inits
    syms = [:θ, :ω, :Pmech]
    values = NetworkDynamics.get_defaults_or_inits(vm, syms)
    @test values[1] == 0.5  # θ has default
    @test values[2] == 0.1  # ω has init
    @test values[3] === nothing  # Pmech has neither

    # Test with Network access
    nw = basenetwork()
    syms = [:θ, :ω, :Pmech]
    # Set up test data
    set_default!(nw, VIndex(1, :θ), 0.5)
    delete_default!(nw, VIndex(1,:ω))
    delete_default!(nw, VIndex(1,:Pmech))
    set_init!(nw, VIndex(1, :ω), 0.1)
    # Test with alternate missing_val
    defaults = NetworkDynamics.get_defaults(nw, VIndex(1), syms; missing_val=:MISSING)
    @test defaults[1] == 0.5  # θ has default
    @test defaults[2] === :MISSING  # ω has no default
    @test defaults[3] === :MISSING  # Pmech has no default
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

@testset "Pattern matching in metadata functions" begin
    # Create a test vertex model
    vm = Lib.swing_mtk()

    # Test happy path - string pattern matching
    set_metadata!(vm, "θ", :test_key, "test_value")
    @test get_metadata(vm, :θ, :test_key) == "test_value"

    # Test happy path - regex pattern matching
    set_default!(vm, r"^ω", 1.5)
    @test get_default(vm, :ω) == 1.5

    # Test no match - should throw ArgumentError
    @test_throws ArgumentError set_metadata!(vm, "nonexistent", :test_key, "value")
    @test_throws ArgumentError set_init!(vm, r"xyz", 2.0)

    set_bounds!(vm, "Pmech", (0.0, 10.0))
    @test get_bounds(vm, :Pmech) == (0.0, 10.0)

    set_guess!(vm, "M", 5.0)
    @test get_guess(vm, :M) == 5.0

    @test_throws ArgumentError set_guess!(vm, "", 5.0)
    @test_throws ArgumentError set_guess!(vm, r"", 5.0)

    # Test basic pattern matching functions
    @test has_metadata(vm, "θ", :test_key)
    @test has_metadata(vm, r"^θ", :test_key)
    @test get_metadata(vm, "θ", :test_key) == "test_value"
    @test get_metadata(vm, r"^θ", :test_key) == "test_value"
    @test delete_metadata!(vm, "θ", :test_key)
    @test !has_metadata(vm, :θ, :test_key)
    @test_throws ArgumentError has_metadata(vm, "nonexistent", :test_key)
    @test_throws ArgumentError get_metadata(vm, r"xyz", :test_key)
end

@testset "set_defaults! function" begin
    nw = basenetwork()
    modified_state = NWState(nw)

    # Modify state values using variable_symbols for correct iteration order
    for (i, sni) in enumerate(variable_symbols(nw))
        if i % 3 == 0
            modified_state[sni] = 2.5 + i/10
        end
    end

    # Modify parameter values using parameter_symbols for correct iteration order
    for (i, sni) in enumerate(parameter_symbols(nw))
        if i % 2 == 0
            modified_state.p[sni] = 3.7 + i/10
        end
    end

    # Test set_defaults! with NWState
    set_defaults!(nw, modified_state)

    # Check if the defaults were updated correctly
    # Need to create a new state to ensure we're reading from the network's defaults
    new_state = NWState(nw)

    for (i, sni) in enumerate(variable_symbols(nw))
        if i % 3 == 0
            @test new_state[sni] ≈ 2.5 + i/10
        end
    end

    for (i, sni) in enumerate(parameter_symbols(nw))
        if i % 2 == 0
            @test new_state.p[sni] ≈ 3.7 + i/10
        end
    end

    # Now test NWParameter separately
    nw2 = basenetwork()
    modified_params = NWParameter(nw2)

    # Modify parameter values using parameter_symbols
    for (i, sni) in enumerate(parameter_symbols(nw2))
        if i % 3 == 1
            modified_params[sni] = 4.2 + i/10
        end
    end

    # Test set_defaults! with NWParameter
    set_defaults!(nw2, modified_params)

    # Check if the defaults were updated correctly
    new_params = NWParameter(nw2)

    for (i, sni) in enumerate(parameter_symbols(nw2))
        if i % 3 == 1
            @test new_params[sni] ≈ 4.2 + i/10
        end
    end

    # Test that NaN values are ignored
    nw3 = basenetwork()
    nan_state = NWState(nw3)

    # Set some values to NaN
    for (i, sni) in enumerate(variable_symbols(nw3))
        if i % 5 == 0
            nan_state[sni] = NaN
        elseif i % 5 == 1
            nan_state[sni] = 7.1
        end
    end

    set_defaults!(nw3, nan_state)

    # Check that NaN values were not set but others were
    new_state3 = NWState(nw3)
    for (i, sni) in enumerate(variable_symbols(nw3))
        if i % 5 == 0
            # The NaN value should not be set, so it should retain original value
            @test nan_state[sni] != new_state3[sni] || isnan(nan_state[sni])
        elseif i % 5 == 1
            @test new_state3[sni] ≈ 7.1
        end
    end
end

@testset "strip_metadata! function" begin
    # Create a test vertex model
    vm = Lib.swing_mtk()

    # Set up various metadata types
    set_default!(vm, :θ, 1.0)
    set_default!(vm, :ω, 2.0)
    set_guess!(vm, :M, 5.0)
    set_init!(vm, :Pmech, 3.0)
    set_metadata!(vm, :θ, :custom, "test_value")

    # Verify metadata was set
    @test has_default(vm, :θ)
    @test has_default(vm, :ω)
    @test has_guess(vm, :M)
    @test has_init(vm, :Pmech)
    @test has_metadata(vm, :θ, :custom)

    # Test stripping defaults
    strip_metadata!(vm, :default)
    @test !has_default(vm, :θ)
    @test !has_default(vm, :ω)
    # Other metadata should remain
    @test has_guess(vm, :M)
    @test has_init(vm, :Pmech)
    @test has_metadata(vm, :θ, :custom)

    # Test stripping custom metadata
    strip_metadata!(vm, :custom)
    @test !has_metadata(vm, :θ, :custom)
    # Other metadata should remain
    @test has_guess(vm, :M)
    @test has_init(vm, :Pmech)
end
