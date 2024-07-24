module ModelingToolkitExt

using ModelingToolkit: Symbolic, istree, operation, arguments, build_function
using ModelingToolkit: Equation, ODESystem, Differential
using ModelingToolkit: full_equations, get_variables, structural_simplify, getname, unwrap
using ModelingToolkit: full_parameters, unknowns,independent_variable, observed, defaults
using ArgCheck: @argcheck
using LinearAlgebra: Diagonal, I

import NetworkDynamics: ODEVertex, StaticEdge

include("ModelingToolkitUtils.jl")

function ODEVertex(sys::ODESystem, inputs, outputs; verbose=false)
    gen = generate_io_function(sys, inputs, outputs; type=:ode, verbose)

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
    ODEVertex(;f=gen.f_ip, sym, def, psym, pdef, depth)
end

function StaticEdge(sys::ODESystem, inputs, outputs; verbose=false)
    gen = generate_io_function(sys, inputs, outputs; type=:static, verbose)

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
    StaticEdge(;f=gen.f_ip, sym, def, psym, pdef, depth)
end

function generate_io_function(_sys::ODESystem, inputs, outputs;
                              expression=Val{false}, verbose=false, type=:auto)

    # f_* may be given in namepsace version or as symbols
    inputs  = _resolve_var.(Ref(_sys), inputs)
    outputs = _resolve_var.(Ref(_sys), outputs)

    sys, _ = structural_simplify(_sys, (inputs, outputs))
    full_equations(sys)
    observed(sys)
    observed_eqs

    @argcheck inputs ⊆ Set(full_parameters(sys))
    @argcheck outputs ⊆ Set(unknowns(sys))

    eqs = full_equations(sys)

    @argcheck isempty(rhs_differentials(eqs)) "RHS should not contain any differentials at this point."

    # make sure that outputs appear first
    states = vcat(outputs, unknowns(sys)) |> unique
    params = setdiff(full_parameters(sys), Set(inputs))

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
        mass_matrix = generate_massmatrix(eqs)
        verbose && @info "Reordered by states and generated mass matrix" mass_matrix
    elseif type == :static
        all_static = all(isequal(:explicit_algebraic), first.(eq_type.(eqs)))
        all_static || throw(ArgumentError("Equations of system are not static!"))
        mass_matrix = nothing
    else
        throw(ArgumentError("Unknown type $type"))
    end

    # now generate the actual functions

    iv = independent_variable(sys)
    formulas = [eq.rhs for eq in eqs]

    if type == :ode
        f_oop, f_ip = build_function(formulas, states, inputs, params, iv; expression)
    elseif type == :static
        f_oop, f_ip = build_function(formulas, inputs, params, iv; expression)
    end

    if !isempty(observed(sys))
        @warn "Function generation for observed states not implemented yet."
    end

    return (;f_oop=f_oop, f_ip=f_ip,
            mass_matrix,
            states,
            inputs,
            params)
end

end
