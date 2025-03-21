module MTKExt

using ModelingToolkit: Symbolic, iscall, operation, arguments, build_function
using ModelingToolkit: ModelingToolkit, Equation, ODESystem, Differential
using ModelingToolkit: full_equations, get_variables, structural_simplify, getname, unwrap
using ModelingToolkit: parameters, unknowns, independent_variables, observed, defaults
using Symbolics: Symbolics, fixpoint_sub, substitute
using RecursiveArrayTools: RecursiveArrayTools
using ArgCheck: @argcheck
using LinearAlgebra: Diagonal, I

using NetworkDynamics: NetworkDynamics, set_metadata!,
                       PureFeedForward, FeedForward, NoFeedForward, PureStateMap
import NetworkDynamics: VertexModel, EdgeModel, AnnotatedSym

include("MTKUtils.jl")

"""
    VertexModel(sys::ODESystem, inputs, outputs;
                verbose=false, name=getname(sys), extin=nothing, ff_to_constraint=true, kwargs...)

Create a `VertexModel` object from a given `ODESystem` created with ModelingToolkit.
You need to provide 2 lists of symbolic names (`Symbol` or `Vector{Symbols}`):
- `inputs`: names of variables in you equation representing the aggregated edge states
- `outputs`: names of variables in you equation representing the node output

Additional kw arguments:
- `name`: Set name of the component model. Will be lifted from the ODESystem name.
- `extin=nothing`: Provide external inputs as pairs, i.e. `extin=[:extvar => VIndex(1, :a)]`
   will bound the variable `extvar(t)` in the equations to the state `a` of the first vertex.
- `ff_to_constraint=true`: Controlls, whether output transformations `g` which depend on inputs should be
  transformed into constraints. Defaults to true since ND.jl does not handle vertices with FF yet.
"""
function VertexModel(sys::ODESystem, inputs, outputs; verbose=false, name=getname(sys),
                     ff_to_constraint=true, extin=nothing, kwargs...)
    warn_events(sys)
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
    EdgeModel(sys::ODESystem, srcin, dstin, AntiSymmetric(dstout); kwargs...)

Create a `EdgeModel` object from a given `ODESystem` created with ModelingToolkit for **single sided models**.

Here you only need to provide one list of output symbols: `dstout`.
To make it clear how to handle the single-sided output definiton, you musst wrap
the symbol vector in
- `AntiSymmetric(dstout)`,
- `Symmetric(dstout)`, or
- `Directed(dstout)`.

Additional `kwargs` are the same as for the double-sided EdgeModel MTK constructor.
"""
EdgeModel(sys::ODESystem, srcin, dstin, dstout; kwargs...) = EdgeModel(sys, srcin, dstin, nothing, dstout; kwargs...)

"""
    EdgeModel(sys::ODESystem, srcin, srcout, dstin, dstout;
              verbose=false, name=getname(sys), extin=nothing, ff_to_constraint=false, kwargs...)

Create a `EdgeModel` object from a given `ODESystem` created with ModelingToolkit.
You need to provide 4 lists of symbolic names (`Symbol` or `Vector{Symbols}`):
- `srcin`: names of variables in you equation representing the node state at the source
- `dstin`: names of variables in you equation representing the node state at the destination
- `srcout`: names of variables in you equation representing the output at the source
- `dstout`: names of variables in you equation representing the output at the destination

Additional kw arguments:
- `name`: Set name of the component model. Will be lifted from the ODESystem name.
- `extin=nothing`: Provide external inputs as pairs, i.e. `extin=[:extvar => VIndex(1, :a)]`
   will bound the variable `extvar(t)` in the equations to the state `a` of the first vertex.
- `ff_to_constraint=false`: Controlls, whether output transformations `g` which depend on inputs should be
  transformed into constraints.
"""
function EdgeModel(sys::ODESystem, srcin, dstin, srcout, dstout; verbose=false, name=getname(sys),
                   ff_to_constraint=false, extin=nothing, kwargs...)
    warn_events(sys)
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
    implicit_outputs = Set{Symbolic}() # fully implicit outputs which do not appear in the equations
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
        _definedoutputs = alloutputs ∩ all_eq_vars
        if !(Set(_definedoutputs) == Set(alloutputs))
            implicit_outputs = setdiff(alloutputs, _definedoutputs)
            verbose && @warn "The specified outputs $implicit_outputs do not appear in the equations of the system!"
        end
        structural_simplify(_sys, (_openinputs, _definedoutputs); simplify=false)[1]
    end

    allparams = parameters(sys) # contains inputs!
    @argcheck allinputs ⊆ Set(allparams) ∪ missing_inputs
    params = setdiff(allparams, Set(allinputs))

    # extract the main equations and observed equations
    eqs::Vector{Equation} = ModelingToolkit.subs_constants(full_equations(sys))
    fix_metadata!(eqs, sys);

    # assert the ordering of states and equations
    explicit_states = Symbolic[eq_type(eq)[2] for eq in eqs if !isnothing(eq_type(eq)[2])]
    implicit_states = setdiff(unknowns(sys), explicit_states) ∪ implicit_outputs
    states = map(eqs) do eq
        type = eq_type(eq)
        isnothing(type[2]) ? pop!(implicit_states) : type[2]
    end

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

    # extract observed equations. They might depend on eachother so resolve them
    obs_subs = Dict(eq.lhs => eq.rhs for eq in observed(sys))
    obseqs = map(observed(sys)) do eq
        expanded_rhs = fixpoint_sub(eq.rhs, obs_subs)
        eq.lhs ~ ModelingToolkit.subs_constants(expanded_rhs)
    end
    fix_metadata!(obseqs, sys);
    # obs can only depend on parameters (including allinputs) or states
    obs_deps = _all_rhs_symbols(obseqs)
    if !(obs_deps ⊆ Set(allparams) ∪ Set(states) ∪ independent_variables(sys))
        @warn "obs_deps !⊆ parameters ∪ unknowns. Difference: $(setdiff(obs_deps, Set(allparams) ∪ Set(states)))"
    end

    # if some states shadow outputs (out ~ state in observed)
    # switch their names. I.e. prioritize use of name `out`
    renamings = Dict()
    for eq in obseqs
        if eq.lhs ∈ Set(alloutputs) && iscall(eq.rhs) &&
            operation(eq.rhs) isa Symbolics.BasicSymbolic && eq.rhs ∈ Set(states)
            verbose && @info "Encountered trivial equation $eq. Swap out $(eq.lhs) <=> $(eq.rhs) everywhere."
            renamings[eq.lhs] = eq.rhs
            renamings[eq.rhs] = eq.lhs
        end
    end
    if !isempty(renamings)
        eqs = map(eq -> substitute(eq, renamings), eqs)
        obseqs = map(eq -> substitute(eq, renamings), obseqs)
        states = map(s -> substitute(s, renamings), states)
        verbose && @info "New States:" states
    end

    # find the output equations, this might remove them from obseqs!
    outeqs = Equation[]
    for out in Iterators.flatten(outputss)
        if out ∈ Set(states)
            push!(outeqs, out ~ out)
        else
            idx = findfirst(eq -> isequal(eq.lhs, out), obseqs)
            if isnothing(idx)
                throw(ArgumentError("Output $out was neither foundin states nor in observed equations."))
            end
            eq = obseqs[idx]
            if !isempty(rhs_differentials(eq))
                println(obs_subs[out])
                throw(ArgumentError("Algebraic FF equation for output $out contains differentials in the RHS: $(rhs_differentials(eq))"))
            end
            deleteat!(obseqs, idx)

            if ff_to_constraint && !isempty(get_variables(eq.rhs) ∩ allinputs)
                verbose && @info "Output $out would lead to FF in g, promote to state instead."
                push!(eqs, 0 ~ eq.lhs - eq.rhs)
                push!(states, eq.lhs)
                push!(outeqs, eq.lhs ~ eq.lhs)
            else
                push!(outeqs, eq)
            end
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

        # create massmatrix, we don't use the method provided by ODESystem because of reordering
        mm = generate_massmatrix(eqs)
        verbose && @info "Generated mass matrix" mm
        mm
    end

    iv = only(independent_variables(sys))
    out_deps = _all_rhs_symbols(outeqs)
    fftype = _determine_fftype(out_deps, states, allinputs, params, iv)

    # filter out unnecessary parameters
    var_deps = _all_rhs_symbols(eqs)
    unused_params = Set(setdiff(params, (var_deps ∪ out_deps))) # do not exclud obs_deps
    if verbose && !isempty(unused_params)
        @info "Parameters $(unused_params) do not appear in equations of f and g and will be marked as unused."
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

    return (;
        f=f_ip, g=g_ip,
        mass_matrix,
        states,
        inputss,
        outputss,
        obsstates,
        fftype,
        obsf = obsf_ip,
        equations=formulas,
        outputeqs=Dict(Iterators.flatten(outputss) .=> gformulas),
        observed=Dict(getname.(obsstates) .=> obsformulas),
        odesystem=sys,
        params,
        unused_params
    )
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
