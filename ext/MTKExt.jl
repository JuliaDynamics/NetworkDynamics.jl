module MTKExt

using ModelingToolkit: Symbolic, iscall, operation, arguments, build_function
using ModelingToolkit: ModelingToolkit, Equation, ODESystem, Differential
using ModelingToolkit: full_equations, get_variables, structural_simplify, getname, unwrap
using ModelingToolkit: full_parameters, unknowns, independent_variables, observed, defaults
using ModelingToolkit.Symbolics: Symbolics, fixpoint_sub
using RecursiveArrayTools: RecursiveArrayTools
using ArgCheck: @argcheck
using LinearAlgebra: Diagonal, I

using NetworkDynamics: NetworkDynamics, set_metadata!,
                       PureFeedForward, FeedForward, NoFeedForward, PureStateMap
import NetworkDynamics: VertexFunction, EdgeFunction, AnnotatedSym

include("MTKUtils.jl")

function VertexFunction(sys::ODESystem, inputs, outputs; verbose=false, name=getname(sys), kwargs...)
    warn_events(sys)
    inputs = inputs isa AbstractVector ? inputs : [inputs]
    outputs = outputs isa AbstractVector ? outputs : [outputs]
    gen = generate_io_function(sys, (inputs,), (outputs,); verbose)

    f = gen.f
    g = gen.g
    obsf = gen.obsf

    _sym = getname.(gen.states)
    sym = [s => _get_metadata(sys, s) for s in _sym]

    _psym = getname.(gen.params)
    psym = [s => _get_metadata(sys, s) for s in _psym]

    _obssym = getname.(gen.obsstates)
    obssym = [s => _get_metadata(sys, s) for s in _obssym]

    _insym = getname.(inputs)
    insym = [s => _get_metadata(sys, s) for s in _insym]

    _outsym = getname.(outputs)
    outsym = [s => _get_metadata(sys, s) for s in _outsym]

    mass_matrix = gen.mass_matrix
    c = VertexFunction(;f, g, sym, insym, outsym, psym, obssym,
            obsf, mass_matrix, ff=gen.fftype, name, allow_output_sym_clash=true, kwargs...)
    set_metadata!(c, :observed, gen.observed)
    set_metadata!(c, :equations, gen.equations)
    set_metadata!(c, :outputeqs, gen.outputeqs)
    c
end

EdgeFunction(sys::ODESystem, srcin, dstin, dstout; kwargs...) = EdgeFunction(sys, srcin, dstin, nothing, dstout; kwargs...)
function EdgeFunction(sys::ODESystem, srcin, dstin, srcout, dstout; verbose=false, name=getname(sys), kwargs...)
    warn_events(sys)
    srcin = srcin isa AbstractVector ? srcin : [srcin]
    dstin = dstin isa AbstractVector ? dstin : [dstin]

    singlesided = isnothing(srcout)
    if singlesided && !(dstout isa AnnotatedSym)
        throw(ArgumentError("If you only provide one output (single sided \
        model), it musst be wrapped either in `AntiSymmetric`, `Symmetric` or \
        `Directed`!"))
    end

    if singlesided
        gwrap = NetworkDynamics.wrapper(dstout)
        dstout = NetworkDynamics.sym(dstout)
        outs = (dstout, )
    else
        srcout = srcout isa AbstractVector ? srcout : [srcout]
        dstout = dstout isa AbstractVector ? dstout : [dstin]
        outs = (srcout, dstout)
    end

    gen = generate_io_function(sys, (srcin, dstin), outs; verbose)

    f = gen.f
    g = singlesided ? gwrap(gen.g; ff=gen.fftype) : gen.g
    obsf = gen.obsf

    _sym = getname.(gen.states)
    sym = [s => _get_metadata(sys, s) for s in _sym]

    _psym = getname.(gen.params)
    psym = [s => _get_metadata(sys, s) for s in _psym]

    _obssym = getname.(gen.obsstates)
    obssym = [s => _get_metadata(sys, s) for s in _obssym]

    _insym_src = getname.(srcin)
    insym_src = [s => _get_metadata(sys, s)  for s in _insym_src]
    _insym_dst = getname.(dstin)
    insym_dst = [s => _get_metadata(sys, s)  for s in _insym_dst]
    insym = (;src=insym_src, dst=insym_dst)

    if singlesided
        _outsym_dst = getname.(dstout)
        outsym = [s => _get_metadata(sys, s)  for s in _outsym_dst]
    else
        _outsym_src = getname.(srcout)
        outsym_src = [s => _get_metadata(sys, s)  for s in _outsym_src]
        _outsym_dst = getname.(dstout)
        outsym_dst = [s => _get_metadata(sys, s)  for s in _outsym_dst]
        outsym = (;src=outsym_src, dst=outsym_dst)
    end

    mass_matrix = gen.mass_matrix
    c = EdgeFunction(;f, g, sym, insym, outsym, psym, obssym,
            obsf, mass_matrix, ff=gen.fftype, name,  allow_output_sym_clash=true, kwargs...)
    set_metadata!(c, :observed, gen.observed)
    set_metadata!(c, :equations, gen.equations)
    set_metadata!(c, :outputeqs, gen.outputeqs)
    c
end

"""
For a given system and name, extract all the relevant meta we want to keep for the component function.
"""
function _get_metadata(sys, name)
    nt = (;)
    sym = try
        getproperty_symbolic(sys, name)
    catch e
        if !endswith(string(name), "ˍt") # known for "internal" derivatives
            @warn "Could not extract metadata for $name $(e.msg)"
        end
        return nt
    end
    alldefaults = defaults(sys)
    if haskey(alldefaults, sym)
        def = alldefaults[sym]
        if def isa Symbolic
            def = fixpoint_sub(def, alldefaults)
        end
        def isa Symbolic && error("Could not resolve default $(ModelingToolkit.getdefault(sym)) for $name")
        nt = (; nt..., default=def)
    end

    # check for guess both in symbol metadata and in guesses of system
    # fixes https://github.com/SciML/ModelingToolkit.jl/issues/3075
    if ModelingToolkit.hasguess(sym) || haskey(ModelingToolkit.guesses(sys), sym)
        guess = if ModelingToolkit.hasguess(sym)
            ModelingToolkit.getguess(sym)
        else
            ModelingToolkit.guesses(sys)[sym]
        end
        if guess isa Symbolic
            guess = fixpoint_sub(def, merge(defaults(sys), guesses(sys)))
        end
        guess isa Symbolic && error("Could not resolve guess $(ModelingToolkit.getguess(sym)) for $name")
        nt = (; nt..., guess=guess)
    end

    if ModelingToolkit.hasbounds(sym)
        nt = (; nt..., bounds=ModelingToolkit.getbounds(sym))
    end
    if ModelingToolkit.hasdescription(sym)
        nt = (; nt..., description=ModelingToolkit.getdescription(sym))
    end
    nt
end

function generate_io_function(_sys, inputss::Tuple, outputss::Tuple;
                              expression=Val{false}, verbose=false)
    # TODO: scalarize vector symbolics/equations?

    # f_* may be given in namepsace version or as symbols, resolve to unnamespaced Symbolic
    inputss = map(inputss) do in
        getproperty_symbolic.(Ref(_sys), in)
    end
    allinputs = reduce(union, inputss)
    outputss = map(outputss) do out
        getproperty_symbolic.(Ref(_sys), out)
    end
    alloutputs = reduce(union, outputss)

    sys = if ModelingToolkit.iscomplete(_sys)
        deepcopy(_sys)
    else
        _openinputs = setdiff(allinputs, Set(full_parameters(_sys)))
        structural_simplify(_sys, (_openinputs, alloutputs); simplify=true)[1]
    end

    states = unknowns(sys)
    allparams = full_parameters(sys) # contains inputs!
    @argcheck allinputs ⊆ Set(allparams)
    params = setdiff(allparams, Set(allinputs))

    # extract the main equations and observed equations
    eqs::Vector{Equation} = full_equations(sys)
    fix_metadata!(eqs, sys);
    # check hat there are no rhs differentials in the equations
    if !isempty(rhs_differentials(eqs))
        diffs = rhs_differentials(eqs)
        buf = IOBuffer()
        println(buf, "Equations contain differentials in their rhs: ", diffs)
        # for (i, eqs) in enumerate(eqs)
        #     if !isempty(rhs_differentials(eqs))
        #         println(buf, " - $i: $eqs")
        #     end
        # end
        throw(ArgumentError(String(take!(buf))))
    end
    # generate mass matrix (this might change the equations)
    mass_matrix = begin
        # equations of form o = f(...) have to be transformed to 0 = f(...) - o
        for (i, eq) in enumerate(eqs)
           if eq_type(eq)[1] == :explicit_algebraic
               eqs[i] = 0 ~ eq.rhs - eq.lhs
           end
        end
        verbose && @info "Transformed algebraic eqs" eqs

        # create massmatrix, we don't use the method provided by ODESystem because of reordering
        mm = generate_massmatrix(eqs)
        verbose && @info "Reordered by states and generated mass matrix" mm
        mm
    end

    # extract observed equations. They might depend on eachother so resolve them
    obs_subs = Dict(eq.lhs => eq.rhs for eq in observed(sys))
    obseqs = map(observed(sys)) do eq
        eq.lhs ~ fixpoint_sub(eq.rhs, obs_subs)
    end
    fix_metadata!(obseqs, sys);
    # obs can only depend on parameters (including allinputs) or states
    obs_deps = _all_rhs_symbols(obseqs)
    if !(obs_deps ⊆ Set(allparams) ∪ Set(states) ∪ independent_variables(sys))
        @warn "obs_deps !⊆ parameters ∪ unknowns. Difference: $(setdiff(obs_deps, Set(allparams) ∪ Set(states)))"
    end

    # find the output equations, this might remove the mfrom obseqs!
    outeqs = map(Iterators.flatten(outputss)) do out
        if out ∈ Set(states)
            out ~ out
        else
            idx = findfirst(eq -> isequal(eq.lhs, out), obseqs)
            if isnothing(idx)
                throw(ArgumentError("Output $out was neither foundin states nor in observed equations."))
            end
            eq = obseqs[idx]
            deleteat!(obseqs, idx)
            obseqs
            eq
        end
    end

    iv = only(independent_variables(sys))
    out_deps = _all_rhs_symbols(outeqs)
    fftype = _determine_fftype(out_deps, states, allinputs, params, iv)

    # filter out unnecessary parameters
    var_deps = _all_rhs_symbols(eqs)
    used_params = params ∩ (var_deps ∪ obs_deps ∪ out_deps)
    if Set(params) != Set(used_params)
        verbose && @info "Remove parameters $(collect(setdiff(params, used_params))) which arn't used in the equations."
        params = used_params
    end

    # TODO: explore Symbolcs/SymbolicUtils CSE
    # now generate the actual functions
    formulas = [eq.rhs for eq in eqs]
    if !isempty(formulas)
        _, f_ip = build_function(formulas, states, inputss..., params, iv; cse=true, expression)
    else
        f_ip = nothing
    end

    gformulas = [eq.rhs for eq in outeqs]
    gformargs = if fftype isa PureFeedForward
        (inputss..., params, iv)
    elseif fftype isa FeedForward
        (states, inputss..., params, iv)
    elseif fftype isa NoFeedForward
        (states, params, iv)
    elseif fftype isa PureStateMap
        (states,)
    end
    _, _g_ip = build_function(gformulas, gformargs...; cse=true, expression)
    # for more thatn 1 output wrap funktion
    g_ip = if length(outputss) == 1
        _g_ip
    else
        MultipleOutputWrapper{fftype, length(outputss), typeof(_g_ip)}(_g_ip)
    end

    # and the observed functions
    obsstates = [eq.lhs for eq in obseqs]
    obsformulas = [eq.rhs for eq in obseqs]
    _, obsf_ip = build_function(obsformulas, states, inputss..., params, iv; cse=true, expression)

    return (;f=f_ip, g=g_ip,
            mass_matrix,
            states,
            inputss,
            outputss,
            obsstates,
            fftype,
            obsf = obsf_ip,
            equations=formulas,
            outputeqs=Dict(Iterators.flatten(outputss) .=> gformulas),
            observed=Dict(obsstates .=> obsformulas),
            params)
end

function _determine_fftype(deps, states, allinputs, params, t)
    if isempty(allinputs ∩ deps) # no ff path
        if isempty(params ∩ deps) && !(t ∈ deps) # no p nor t
            return PureStateMap()
        else # needs p or t
            return NoFeedForward()
        end
    else # ff path
        if isempty(states ∩ deps) # without state dependency
            return PureFeedForward()
        else # with state dependency
            return FeedForward()
        end
    end
end

_all_rhs_symbols(eqs) = mapreduce(eq->get_variables(eq.rhs), ∪, eqs, init=Set{Symbolic}())

using PrecompileTools: @setup_workload, @compile_workload
@setup_workload begin
    @compile_workload begin
        # include("precompile_workload.jl")
    end
end

end
