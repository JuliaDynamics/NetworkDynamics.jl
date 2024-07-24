module ModelingToolkitExt

using ModelingToolkit: Symbolic, istree, operation, arguments, build_function
using ModelingToolkit: ModelingToolkit, Equation, ODESystem, Differential
using ModelingToolkit: full_equations, get_variables, structural_simplify, getname, unwrap
using ModelingToolkit: full_parameters, unknowns, independent_variable, observed, defaults
using ModelingToolkit: get_substitutions
using ModelingToolkit.Symbolics: fixpoint_sub
using ArgCheck: @argcheck
using LinearAlgebra: Diagonal, I

using NetworkDynamics: Fiducial
import NetworkDynamics: ODEVertex, StaticEdge

include("ModelingToolkitUtils.jl")

function ODEVertex(sys::ODESystem, inputs, outputs; verbose=false)
    warn_events(sys)
    inputs = inputs isa AbstractVector ? inputs : [inputs]
    outputs = outputs isa AbstractVector ? outputs : [outputs]
    gen = generate_io_function(sys, (inputs, ), outputs; type=:ode, verbose)

    f = gen.f_ip
    defdict = defaults(sys)
    sym = getname.(gen.states)
    def = map(gen.states) do s
        get(defdict, s, nothing)
    end
    psym = getname.(gen.params)
    pdef = map(gen.params) do p
        get(defdict, p, nothing)
    end
    depth = length(outputs)
    obsf = gen.g_ip
    obssym = getname.(gen.params)
    mass_matrix = gen.mass_matrix
    ODEVertex(;f, sym, def, psym, pdef, depth, obssym, obsf, mass_matrix)
end

function StaticEdge(sys::ODESystem, srcin, dstin, outputs, coupling; verbose=false)
    warn_events(sys)
    srcin = srcin isa AbstractVector ? srcin : [srcin]
    dstin = dstin isa AbstractVector ? dstin : [dstin]
    outputs = outputs isa AbstractVector ? outputs : [outputs]
    gen = generate_io_function(sys, (srcin, dstin), outputs; type=:static, verbose)

    f = gen.f_ip
    defdict = defaults(sys)
    sym = getname.(gen.states)
    def = map(gen.states) do s
        get(defdict, s, nothing)
    end
    psym = getname.(gen.params)
    pdef = map(gen.params) do p
        get(defdict, p, nothing)
    end
    obsf = gen.g_ip
    obssym = getname.(gen.params)
    depth = coupling isa Fiducial ? Int(length(outputs)/2) : length(outputs)
    StaticEdge(;f, sym, def, psym, pdef, depth, obssym, obsf, coupling)
end

function generate_io_function(_sys, inputss::Tuple, outputs;
                              expression=Val{false}, verbose=false, type=:auto)
    # TODO: scalarize vector symbolics/equations?

    # f_* may be given in namepsace version or as symbols, resolve to unnamespaced Symbolic
    inputss = map(inputss) do in
        _resolve_var.(Ref(_sys), in)
    end
    allinputs = reduce(union, inputss)
    outputs = _resolve_var.(Ref(_sys), outputs)
    sys, _ = structural_simplify(_sys, (allinputs, outputs); simplify=true)

    # extract the main equations and observed equations
    eqs::Vector{Equation} = full_equations(sys)

    # extract observed equations. They might depend on eachother so resolve them
    obs_subs = Dict(eq.lhs => eq.rhs for eq in observed(sys))
    obseqs = map(observed(sys)) do eq
        eq.lhs ~ fixpoint_sub(eq.rhs, obs_subs)
    end
    # obs can only depend on parameters (including allinputs) or states
    obs_deps = mapreduce(eq -> get_variables(eq.rhs), union, obseqs, init=Symbolic[])
    @assert obs_deps ⊆ Set(full_parameters(sys)) ∪ Set(unknowns(sys))

    @argcheck allinputs ⊆ Set(full_parameters(sys))

    if !(outputs ⊆ Set(unknowns(sys)))
        # structural simplify might remove explicit equations for outputs
        # so we reconstruct equations for them based on the substitutions
        missingouts = setdiff(Set(outputs), unknowns(sys))

        subeqs = get_substitutions(sys).subs
        subs = Dict(eq.lhs => eq.rhs for eq in subeqs)

        @argcheck missingouts ⊆ keys(subs)

        for mout in missingouts
            eq = mout ~ fixpoint_sub(mout, subs)
            push!(eqs, eq)
        end

        # we promoted the missing outputs to "states" again, so we need to remove them from obseqs
        filter!(obseqs) do eq
            eq.lhs ∉ missingouts
        end
    end

    @argcheck isempty(rhs_differentials(eqs)) "RHS should not contain any differentials at this point."

    # make sure that outputs appear first
    states = vcat(outputs, unknowns(sys)) |> unique
    params = setdiff(full_parameters(sys), Set(allinputs))

    eqs = reorder_by_states(eqs, states)

    verbose && @info "Reordered eqs" eqs states

    if type == :auto
        all_static = all(isequal(:explicit_algebraic), first.(eq_type.(eqs)))
        type = all_static ? :static : :ode
        verbose && @info "auto-equation type: $type"
    end

    mass_matrix = if type == :ode
        # equations of form o = f(...) have to be transformed to 0 = f(...) - o
        for (i, eq) in enumerate(eqs)
           if eq_type(eq)[1] == :explicit_algebraic
               eqs[i] = 0 ~ eq.rhs - eq.lhs
           end
        end
        verbose && @info "Transformed algebraic eqs" eqs

        # create massmatrix, we don't use the method provided by ODESystem because of reordering
        mm = generate_massmatrix(eqs)
        verbose && @info "Reordered by states and generated mass matrix" mass_matrix
        mm
    elseif type == :static
        all_static = all(isequal(:explicit_algebraic), first.(eq_type.(eqs)))
        all_static || throw(ArgumentError("Equations of system are not static!"))
        nothing
    else
        throw(ArgumentError("Unknown type $type"))
    end

    # now generate the actual functions
    iv = independent_variable(sys)

    formulas = [eq.rhs for eq in eqs]
    if type == :ode
        f_oop, f_ip = build_function(formulas, states, inputss..., params, iv; expression)
    elseif type == :static
        f_oop, f_ip = build_function(formulas, inputss..., params, iv; expression)
    end

    # and the observed functions
    obsstates = [eq.lhs for eq in obseqs]
    obsformulas = [eq.rhs for eq in obseqs]
    g_oop, g_ip = build_function(obsformulas, states, inputss..., params, iv; expression)

    return (;f_oop, f_ip,
            mass_matrix,
            states,
            inputss,
            obsstates,
            g_oop, g_ip,
            params)
end

end
