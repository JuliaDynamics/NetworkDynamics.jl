module NetworkDynamicsMTKExt

using ModelingToolkit: Symbolic, iscall, operation, arguments, build_function
using ModelingToolkit: ModelingToolkit, Equation, System, Differential
using ModelingToolkit: equations, full_equations, get_variables, mtkcompile, getname, unwrap
using ModelingToolkit: parameters, unknowns, independent_variables, observed, defaults
using Symbolics: Symbolics, fixpoint_sub, substitute
using RecursiveArrayTools: RecursiveArrayTools
using ArgCheck: @argcheck
using LinearAlgebra: Diagonal, I
using SymbolicUtils.Code: Let, Assignment
using SymbolicUtils: SymbolicUtils
using OrderedCollections: OrderedDict

using NetworkDynamics: NetworkDynamics, set_metadata!,
                       PureFeedForward, FeedForward, NoFeedForward, PureStateMap
import NetworkDynamics: VertexModel, EdgeModel, AnnotatedSym

include("MTKExt_utils.jl")

import NetworkDynamics: implicit_output
ModelingToolkit.@register_symbolic implicit_output(x)

"""
    VertexModel(sys::System, inputs, outputs;
                verbose=false, name=getname(sys), extin=nothing, ff_to_constraint=true, kwargs...)

Create a `VertexModel` object from a given `System` created with ModelingToolkit.
You need to provide 2 lists of symbolic names (`Symbol` or `Vector{Symbols}`):
- `inputs`: names of variables in you equation representing the aggregated edge states
- `outputs`: names of variables in you equation representing the node output

Additional kw arguments:
- `name`: Set name of the component model. Will be lifted from the System name.
- `extin=nothing`: Provide external inputs as pairs, i.e. `extin=[:extvar => VIndex(1, :a)]`
   will bound the variable `extvar(t)` in the equations to the state `a` of the first vertex.
- `ff_to_constraint=true`: Controls, whether output transformations `g` which depend on inputs should be
  transformed into constraints. Defaults to true since ND.jl does not handle vertices with FF yet.
"""
function VertexModel(sys::System, inputs, outputs; verbose=false, name=getname(sys),
                     ff_to_constraint=true, extin=nothing, kwargs...)
    warn_missing_features(sys)
    inputs = inputs isa AbstractVector ? inputs : [inputs]
    outputs = outputs isa AbstractVector ? outputs : [outputs]

    if isnothing(extin)
        extin_nwidx = nothing
        ins = (inputs, )
    else
        extin_sym, extin_nwidx = _split_extin(extin)
        ins = (inputs, extin_sym)
    end

    gen = generate_io_function(sys, ins, (outputs,); verbose, ff_to_constraint)

    f = gen.f
    g = gen.g
    obsf = gen.obsf

    _sym = getname.(gen.states)
    sym = [s => _get_metadata(sys, s) for s in _sym]

    _psym = getname.(gen.params)
    psym = map(gen.params) do p
        pname = getname(p)
        md = _get_metadata(sys, pname)
        if p in gen.unused_params
            md = (; md..., unused=true)
        end
        pname => md
    end

    _obssym = getname.(gen.obsstates)
    obssym = [s => _get_metadata(sys, s) for s in _obssym]

    _insym = getname.(inputs)
    insym = [s => _get_metadata(sys, s) for s in _insym]

    _outsym = getname.(outputs)
    outsym = [s => _get_metadata(sys, s) for s in _outsym]

    mass_matrix = gen.mass_matrix
    c = VertexModel(;f, g, sym, insym, outsym, psym, obssym,
            obsf, mass_matrix, ff=gen.fftype, name, extin=extin_nwidx,
            allow_output_sym_clash=true, kwargs...)
    set_metadata!(c, :observed, gen.observed)
    set_metadata!(c, :equations, gen.equations)
    set_metadata!(c, :outputeqs, gen.outputeqs)
    set_metadata!(c, :odesystem, gen.odesystem)
    c
end

"""
    EdgeModel(sys::System, srcin, dstin, AntiSymmetric(dstout); kwargs...)

Create a `EdgeModel` object from a given `System` created with ModelingToolkit for **single sided models**.

Here you only need to provide one list of output symbols: `dstout`.
To make it clear how to handle the single-sided output definition, you must wrap
the symbol vector in
- `AntiSymmetric(dstout)`,
- `Symmetric(dstout)`, or
- `Directed(dstout)`.

Additional `kwargs` are the same as for the double-sided EdgeModel MTK constructor.
"""
EdgeModel(sys::System, srcin, dstin, dstout; kwargs...) = EdgeModel(sys, srcin, dstin, nothing, dstout; kwargs...)

"""
    EdgeModel(sys::System, srcin, dstin, srcout, dstout;
              verbose=false, name=getname(sys), extin=nothing, ff_to_constraint=false, kwargs...)

Create a `EdgeModel` object from a given `System` created with ModelingToolkit.
You need to provide 4 lists of symbolic names (`Symbol` or `Vector{Symbols}`):
- `srcin`: names of variables in you equation representing the node state at the source
- `dstin`: names of variables in you equation representing the node state at the destination
- `srcout`: names of variables in you equation representing the output at the source
- `dstout`: names of variables in you equation representing the output at the destination

Additional kw arguments:
- `name`: Set name of the component model. Will be lifted from the System name.
- `extin=nothing`: Provide external inputs as pairs, i.e. `extin=[:extvar => VIndex(1, :a)]`
   will bound the variable `extvar(t)` in the equations to the state `a` of the first vertex.
- `ff_to_constraint=false`: Controls, whether output transformations `g` which depend on inputs should be
  transformed into constraints.
"""
function EdgeModel(sys::System, srcin, dstin, srcout, dstout; verbose=false, name=getname(sys),
                   ff_to_constraint=false, extin=nothing, kwargs...)
    warn_missing_features(sys)
    srcin = srcin isa AbstractVector ? srcin : [srcin]
    dstin = dstin isa AbstractVector ? dstin : [dstin]

    if isnothing(extin)
        extin_nwidx = nothing
        ins = (srcin, dstin)
    else
        extin_sym, extin_nwidx = _split_extin(extin)
        ins = (srcin, dstin, extin_sym)
    end

    singlesided = isnothing(srcout)
    if singlesided && !(dstout isa AnnotatedSym)
        throw(ArgumentError("If you only provide one output (single sided \
        model), it must be wrapped either in `AntiSymmetric`, `Symmetric` or \
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

    gen = generate_io_function(sys, ins, outs; verbose, ff_to_constraint)

    f = gen.f
    g = singlesided ? gwrap(gen.g) : gen.g
    obsf = gen.obsf

    _sym = getname.(gen.states)
    sym = [s => _get_metadata(sys, s) for s in _sym]

    _psym = getname.(gen.params)
    psym = map(gen.params) do p
        pname = getname(p)
        md = _get_metadata(sys, pname)
        if p in gen.unused_params
            md = (; md..., unused=true)
        end
        pname => md
    end

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
    c = EdgeModel(;f, g, sym, insym, outsym, psym, obssym,
            obsf, mass_matrix, ff=gen.fftype, name, extin=extin_nwidx,
            allow_output_sym_clash=true, kwargs...)
    set_metadata!(c, :observed, gen.observed)
    set_metadata!(c, :equations, gen.equations)
    set_metadata!(c, :outputeqs, gen.outputeqs)
    set_metadata!(c, :odesystem, gen.odesystem)
    c
end

"""
For a given system and name, extract all the relevant meta we want to keep for the component model.
"""
function _get_metadata(sys, name)
    nt = (;)
    sym = try
        getproperty_symbolic(sys, name; might_contain_toplevel_ns=false)
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

        if def isa Symbolic
            # do nothing, as the warning can get annoying
            # @warn "Could not resolve rhs for default term $name = $(ModelingToolkit.getdefault(sym)). Some rhs symbols might not have default values. Leave free."
        elseif def == ModelingToolkit.NoValue || def isa ModelingToolkit.NoValue
            # skip NoValue thing
        else
            nt = (; nt..., default=def)
        end
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

function _split_extin(extin)
    try
        extin_sym   = first.(extin)
        extin_nwidx = last.(extin)
        @assert extin_sym isa Vector{Symbol}
        @assert extin_nwidx isa Vector{<:NetworkDynamics.SymbolicIndex}
        return extin_sym, extin_nwidx
    catch e
        @error "Could not evaluate extin keyword argument!"
    end
end

function generate_io_function(_sys, inputss::Tuple, outputss::Tuple;
                              expression=Val{false}, verbose=false,
                              ff_to_constraint)
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

    missing_inputs = Set{Symbolic}()
    sys = if ModelingToolkit.iscomplete(_sys)
        deepcopy(_sys)
    else
        _openinputs = setdiff(allinputs, Set(parameters(_sys)))
        all_eq_vars = mapreduce(get_variables, union, full_equations(_sys), init=Set{Symbolic}())
        if !(_openinputs ⊆ all_eq_vars)
            missing_inputs = setdiff(_openinputs, all_eq_vars)
            verbose && @warn "The specified inputs ($missing_inputs) do not appear in the equations of the system!"
            _openinputs = setdiff(_openinputs, missing_inputs)
        end

        implicit_outputs = setdiff(alloutputs, all_eq_vars)
        if !isempty(implicit_outputs)
            throw(
                ArgumentError("The outputs $(getname.(implicit_outputs)) do not appear in the equations of the system! \
                    Try to to make them explicit using `implicit_output`\n" *
                    NetworkDynamics.implicit_output_docstring)
            )
        end

        verbose && @info "Simplifying system with inputs $_openinputs and outputs $alloutputs"
        try
            mtkcompile(_sys; inputs=_openinputs, outputs=alloutputs, simplify=false)
        catch e
            if e isa ModelingToolkit.ExtraEquationsSystemException
                msg = "The system could not be compiled because of extra equations! \
                       Sometimes, this can be related to fully implicit output equations. \
                       Check `@doc implicit_output` for more information."
                throw(ArgumentError(msg))
            end
            rethrow(e)
        end
    end

    allparams = parameters(sys) # contains inputs!
    @argcheck allinputs ⊆ Set(allparams) ∪ missing_inputs
    params = setdiff(allparams, Set(allinputs))

    # extract the main equations and observed equations
    eqs::Vector{Equation} = equations(sys)
    obseqs_sorted::Vector{Equation} = observed(sys)
    fix_metadata!(eqs, sys);
    fix_metadata!(obseqs_sorted, sys);

    # get rid of the implicit_output(⋅) terms
    remove_implicit_output_fn!(eqs)
    remove_implicit_output_fn!(obseqs_sorted)

    # assert the ordering of states and equations
    explicit_states = Symbolic[eq_type(eq)[2] for eq in eqs if !isnothing(eq_type(eq)[2])]
    implicit_states = setdiff(unknowns(sys), explicit_states)

    if length(explicit_states) + length(implicit_states) !== length(eqs)
        buf = IOBuffer()
        println(buf, "The number of states does not match the number of equations.")
        println(buf, "Explicit states: ", explicit_states)
        println(buf, "Implicit states: ", implicit_states)
        println(buf, "$(length(eqs)) Equations.")
        throw(ArgumentError(String(take!(buf))))
    end

    states = map(eqs) do eq
        type = eq_type(eq)
        isnothing(type[2]) ? pop!(implicit_states) : type[2]
    end

    # check hat there are no rhs differentials in the equations
    if !isempty(rhs_differentials(vcat(eqs, obseqs_sorted)))
        diffs = rhs_differentials(vcat(eqs, obseqs_sorted))
        buf = IOBuffer()
        println(buf, "Equations contain differentials in their rhs: ", diffs)
        # for (i, eqs) in enumerate(eqs)
        #     if !isempty(rhs_differentials(eqs))
        #         println(buf, " - $i: $eqs")
        #     end
        # end
        throw(ArgumentError(String(take!(buf))))
    end

    # obs can only depend on parameters (including allinputs) or states
    obs_subs = OrderedDict(eq.lhs => eq.rhs for eq in obseqs_sorted)
    obs_deps = _all_rhs_symbols(fixpoint_sub(obseqs_sorted, obs_subs))
    if !(obs_deps ⊆ Set(allparams) ∪ Set(states) ∪ independent_variables(sys))
        @warn "obs_deps !⊆ parameters ∪ unknowns. Difference: $(setdiff(obs_deps, Set(allparams) ∪ Set(states)))"
    end

    # if some outputs are direct aliases for states
    # switch their names. I.e. prioritize use of name `out`
    function state_alias(output)
        if output ∈ keys(obs_subs) &&
           iscall(obs_subs[output]) &&
           operation(obs_subs[output]) isa Symbolics.BasicSymbolic
            return state_alias(obs_subs[output])
        elseif output ∈ Set(states)
            return output
        else
            return nothing
        end
    end

    renamings = Dict()
    for output in alloutputs
        sa = state_alias(output)
        if !isnothing(sa)
            verbose && @info "Output $output ≙ $sa, prioritize output name over state name."
            renamings[output] = sa
            renamings[sa] = output
        end
    end
    if !isempty(renamings)
        eqs = map(eq -> substitute(eq, renamings), eqs)
        obseqs_sorted = map(eq -> substitute(eq, renamings), obseqs_sorted)
        obs_subs = OrderedDict(eq.lhs => eq.rhs for eq in obseqs_sorted)
        states = map(s -> substitute(s, renamings), states)
        verbose && @info "New states with applied output aliases:" states
    end

    # find the output equations, this might remove them from obseqs_sorted (obs_subs stays intact)
    outeqs = Equation[]
    for out in Iterators.flatten(outputss)
        if out ∈ Set(states)
            push!(outeqs, out ~ out)
        elseif out ∈ keys(obs_subs)
            # if its a observed, we need to check for ff behavior
            fulleq = out ~ fixpoint_sub(obs_subs[out], obs_subs)
            if ff_to_constraint && !isempty(get_variables(fulleq.rhs) ∩ allinputs)
                verbose && @info "Output $out would lead to FF in g, promote to state instead."
                # not observed anymore, delete from observed and put in equations
                push!(eqs, 0 ~ out - obs_subs[out])
                push!(states, out)
                push!(outeqs, out ~ out)

                deleteat!(obseqs_sorted, findfirst(eq -> isequal(eq.lhs, out), obseqs_sorted))
                delete!(obs_subs, out)
            else # "normal" observed state
                push!(outeqs, out ~ obs_subs[out])
                # delete from obs equations but *not* from obs_subs (otherwise can't be reference)
                # in equations
                deleteat!(obseqs_sorted, findfirst(eq -> isequal(eq.lhs, out), obseqs_sorted))
            end
        else
            throw(ArgumentError("Output $out was neither foundin states nor in observed equations."))
        end
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

        # create massmatrix, we don't use the method provided by System because of reordering
        mm = generate_massmatrix(eqs)
        verbose && @info "Generated mass matrix" mm
        mm
    end

    iv = only(independent_variables(sys))

    out_deps = _all_rhs_symbols(fixpoint_sub(outeqs, obs_subs))
    fftype = _determine_fftype(out_deps, states, allinputs, params, iv)

    # filter out unnecessary parameters
    unused_params = let
        # we need to collect obs and var deps separatly, because replacenment of obs in vars might lead to
        # symbolic simplifications which we don't have in the actual equations later!
        # i.e. dt(x) ~ x - x0 and x ~ x0 + y leads to dt(x) ~ y in substitution but is not resolved in actual f
        var_deps = _all_rhs_symbols(eqs)
        obs_deps = _all_rhs_symbols(obs_subs)
        Set(setdiff(params, (var_deps ∪ obs_deps ∪ out_deps))) # do not exclud obs_deps
    end
    if verbose && !isempty(unused_params)
        @info "Parameters $(unused_params) do not appear in equations of f and g and will be marked as unused."
    end

    # TODO: explore Symbolcs/SymbolicUtils CSE
    # now generate the actual functions
    if !isempty(eqs)
        formulas = _get_formulas(eqs, obs_subs)
        _, f_ip = build_function(formulas, states, inputss..., params, iv; cse=false, expression)
    else
        f_ip = nothing
    end

    # find all observable assigments necessary for outeqs
    gformulas = _get_formulas(outeqs, obs_subs)
    gformargs = if fftype isa PureFeedForward
        (inputss..., params, iv)
    elseif fftype isa FeedForward
        (states, inputss..., params, iv)
    elseif fftype isa NoFeedForward
        (states, params, iv)
    elseif fftype isa PureStateMap
        (states,)
    end
    _, _g_ip = build_function(gformulas, gformargs...; cse=false, expression)
    # for more thatn 1 output wrap funktion
    g_ip = if length(outputss) == 1
        _g_ip
    else
        MultipleOutputWrapper{fftype, length(outputss), typeof(_g_ip)}(_g_ip)
    end

    obsstates = [eq.lhs for eq in obseqs_sorted]
    if !isempty(obsstates)
        obsformulas = _get_formulas([s ~ s for s in obsstates], obs_subs)
        _, obsf_ip = build_function(obsformulas, states, inputss..., params, iv; cse=false, expression)
    else
        obsf_ip = nothing
    end

    return (;
        f=f_ip, g=g_ip,
        mass_matrix,
        states,
        inputss,
        outputss,
        obsstates,
        fftype,
        obsf = obsf_ip,
        equations=eqs,
        outputeqs=outeqs,
        observed=obseqs_sorted,
        odesystem=sys,
        params,
        unused_params
    )
end

function _get_formulas(eqs, obs_subs)
    # Bit hacky, were building a function like this,
    # where all (necessary) obs and eqs are contained in the bgin block of the first output
    # out[1] = begin
    #     obs1   = ...
    #     obs2   = ...
    #     ...
    #     state1 = ...
    #     state2 = ...
    #     ...
    #     state1 # ens up in out[1]
    # end
    # out[2] = state2
    # ...
    isempty(eqs) && return []
    obsdeps = _collect_deps_on_obs([eq.rhs for eq in eqs], obs_subs)
    obs_assignments = [Assignment(k, v) for (k,v) in obs_subs if k ∈ obsdeps]

    # implicit equations are not use via assigments, so we filter for e
    eqs_assignments = [Assignment(eq.lhs, eq.rhs) for eq in eqs
                          if !isequal(eq.lhs, eq.rhs) && !isequal(eq.lhs, 0)]
    # since implicit eqs did not end up in assighmets, we use the rhs
    out = [isequal(eq.lhs, 0) ? eq.rhs : eq.lhs for eq in eqs]

    [Let(vcat(obs_assignments, eqs_assignments), out[1], false), out[2:end]...]
end
function _collect_deps_on_obs(terms, obs_subs)
    deps = Set{Symbolic}()
    for term in terms
        _collect_deps_on_obs!(deps, obs_subs, term)
    end
    deps
end
function _collect_deps_on_obs!(deps, obs_subs, term)
    termdeps = get_variables(term)
    for sym in termdeps
        if haskey(obs_subs, sym)
            # check recursively whether the observed depends on other observed
            _collect_deps_on_obs!(deps, obs_subs, obs_subs[sym])
            push!(deps, sym)
        end
    end
    deps
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

_all_rhs_symbols(eqs) = mapreduce(eq->get_variables(eq isa Pair ? eq.second : eq.rhs), ∪, eqs, init=Set{Symbolic}())

using PrecompileTools: @compile_workload
@compile_workload begin
    # include("MTKExt_precomp_workload.jl")
end

end
