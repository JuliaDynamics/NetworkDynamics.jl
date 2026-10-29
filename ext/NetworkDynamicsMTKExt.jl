module NetworkDynamicsMTKExt

using ModelingToolkitBase: iscall, operation, arguments, build_function
using ModelingToolkitBase: ModelingToolkitBase, Equation, System, Differential
using ModelingToolkitBase: equations, full_equations, get_variables, mtkcompile, getname, unwrap
using ModelingToolkitBase: parameters, unknowns, independent_variables, observed, initial_conditions
using Symbolics: Symbolics, fixpoint_sub, substitute, LinearExpander
using RecursiveArrayTools: RecursiveArrayTools
using ArgCheck: @argcheck
using Graphs: Graphs
using Hungarian: Hungarian, hungarian
using LinearAlgebra: LinearAlgebra, Diagonal, I, UniformScaling
using SymbolicUtils.Code: Let, Assignment, Func, DestructuredArgs, unwrap_const
using Moshi: Moshi
using Moshi.Match: @match
using Accessors: @reset
using NetworkDynamics: MacroTools
using SymbolicUtils: SymbolicUtils
using OrderedCollections: OrderedDict
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII

using NetworkDynamics: NetworkDynamics, set_metadata!, ComponentPostprocessing,
                       PureFeedForward, FeedForward, NoFeedForward, PureStateMap,
                       MultipleOutputWrapper, inline_repr, multiline_repr,
                       AliasMap,
                       ResolutionRule, set_obsrules!,
                       assert_initformula_compat, assert_guessformula_compat
import NetworkDynamics: VertexModel, EdgeModel, AnnotatedSym, InitFormula, add_initformula!, GuessFormula, add_guessformula!

const ST = SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}

include("MTKExt_utils.jl")
include("MTKExt_codegen.jl")
include("MTKExt_simplification.jl")
include("MTKExt_cache.jl")
include("MTKExt_cached_component.jl")

import NetworkDynamics: implicit_output, RHSDifferentialsError
Symbolics.@register_symbolic implicit_output(x)

# `System`s and `Equation`s are large symbolic compile artifacts that are treated as
# immutable once stored in component metadata (e.g. `:odesystem`, `:equations`). Opt
# them into sharing so copying/reconstructing a component does not deep-copy them (which
# is orders of magnitude slower). See `NetworkDynamics.dealias_metadata`.
NetworkDynamics.dealias_metadata(x::Union{System,Equation}) = x

const MTKCOMPILE_DEFAULT = Ref{Union{Symbol,Bool}}(false)
function NetworkDynamics.set_mtkcompile!(value)
    MTKCOMPILE_DEFAULT[] = value
end

"""
    VertexModel(sys::System, inputs, outputs;
                verbose=false, name=getname(sys), extin=nothing, ff_to_constraint=true,
                mtkcompile=nothing, kwargs...)

Create a `VertexModel` object from a given `System` created with ModelingToolkitBase.
You need to provide 2 lists of symbolic names (`Symbol` or `Vector{Symbols}`):
- `inputs`: names of variables in your equation representing the aggregated edge states
- `outputs`: names of variables in your equation representing the node output

Additional kw arguments:
- `name`: Set name of the component model. Will be lifted from the System name.
- `extin=nothing`: Provide external inputs as pairs, i.e. `extin=[:extvar => VIndex(1, :a)]`
   will bound the variable `extvar(t)` in the equations to the state `a` of the first vertex.
- `ff_to_constraint=true`: Controls, whether output transformations `g` which depend on inputs should be
  transformed into constraints. Defaults to true since ND.jl does not handle vertices with FF yet.
- `assume_io_coupling=false`: If true, assumes direct dependency chain from outputs to inputs by adding
  implicit output terms to all inputs. This forces MTK to consider the dependency during compilation
  and can help in cases where MTK simplification results in derivatives of inputs.
- `mtkcompile=nothing`: Controls the simplification backend. `false` (default) uses NetworkDynamics'
  own algebraic reduction pipeline. `true` delegates to MTK's `mtkcompile`: when `ModelingToolkit` is
  loaded this performs full structural simplification; when only `ModelingToolkitBase` is available it
  falls back to the stub `mtkcompile` (no structural simplification). `nothing` inherits the global
  default set via `NetworkDynamics.set_mtkcompile!`. `:compare` runs both and prints a comparison
  (for debugging).
"""
function VertexModel(
    sys::System, inputs, outputs;
    verbose=false,
    name=getname(sys),
    ff_to_constraint=true,
    assume_io_coupling=false,
    extin=nothing,
    mtkcompile=nothing,
    kwargs...
)
    mtkcompile = isnothing(mtkcompile) ? MTKCOMPILE_DEFAULT[] : mtkcompile
    if mtkcompile==:compare
        args = (sys, inputs, outputs)
        kwargs = (; verbose, name, ff_to_constraint, assume_io_coupling, extin, kwargs...)
        return _compare_mtkcompile(VertexModel, args, kwargs)
    end

    warn_missing_features(sys)
    sys = resolve_bound_to(sys)
    inputs = inputs isa AbstractVector ? inputs : [inputs]
    outputs = outputs isa AbstractVector ? outputs : [outputs]

    if isnothing(extin)
        extin_nwidx = nothing
        ins = (inputs, )
    else
        extin_sym, extin_nwidx = _split_extin(extin)
        ins = (inputs, extin_sym)
    end

    gen = generate_io_function_cached(sys, ins, (outputs,); verbose, ff_to_constraint, assume_io_coupling, mtkcompile)

    f = gen.f
    g = gen.g
    obsf = gen.obsf

    sysdefaults = initial_conditions(sys)
    sysguesses = ModelingToolkitBase.guesses(sys)
    _sym = getname.(gen.states)
    sym = [s => _get_metadata(sys, s, sysdefaults, sysguesses) for s in _sym]

    _psym = getname.(gen.params)
    psym = map(gen.params) do p
        pname = getname(p)
        md = _get_metadata(sys, pname, sysdefaults, sysguesses)
        if p in gen.unused_params
            md[:unused] = true
        end
        pname => md
    end

    _obssym = getname.(gen.obsstates)
    obssym = [s => _get_metadata(sys, s, sysdefaults, sysguesses) for s in _obssym]

    _insym = getname.(inputs)
    insym = [s => _get_metadata(sys, s, sysdefaults, sysguesses) for s in _insym]

    _outsym = getname.(outputs)
    outsym = [s => _get_metadata(sys, s, sysdefaults, sysguesses) for s in _outsym]

    mass_matrix = gen.mass_matrix
    c = VertexModel(;f, g, sym, insym, outsym, psym, obssym,
            obsf, mass_matrix, ff=gen.fftype, name, extin=extin_nwidx,
            aliasmap=gen.aliasmap, allow_output_sym_clash=true, kwargs...)
    set_metadata!(c, :observed, gen.observed)
    set_obsrules!(c, gen.obsrules)
    set_metadata!(c, :equations, gen.equations)
    set_metadata!(c, :full_equations, gen.full_equations)
    set_metadata!(c, :outputeqs, gen.outputeqs)
    set_metadata!(c, :full_outputeqs, gen.full_outputeqs)
    set_metadata!(c, :odesystem, sys)

    for formula in something(gen.initformulas, ())
        add_initformula_lenient!(c, formula)
    end
    for formula in something(gen.guessformulas, ())
        add_guessformula_lenient!(c, formula)
    end

    # apply postprocessing functions from subcomponent metadata
    apply_component_postprocessing!(c)
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
              verbose=false, name=getname(sys), extin=nothing, ff_to_constraint=false,
              mtkcompile=nothing, kwargs...)

Create a `EdgeModel` object from a given `System` created with ModelingToolkitBase.
You need to provide 4 lists of symbolic names (`Symbol` or `Vector{Symbols}`):
- `srcin`: names of variables in your equation representing the node state at the source
- `dstin`: names of variables in your equation representing the node state at the destination
- `srcout`: names of variables in your equation representing the output at the source
- `dstout`: names of variables in your equation representing the output at the destination

Additional kw arguments:
- `name`: Set name of the component model. Will be lifted from the System name.
- `extin=nothing`: Provide external inputs as pairs, i.e. `extin=[:extvar => VIndex(1, :a)]`
   will bound the variable `extvar(t)` in the equations to the state `a` of the first vertex.
- `ff_to_constraint=false`: Controls, whether output transformations `g` which depend on inputs should be
  transformed into constraints.
- `assume_io_coupling=false`: If true, assumes direct dependency chain from outputs to inputs by adding
  implicit output terms to all inputs. This forces MTK to consider the dependency during compilation
  and can help in cases where MTK simplification results in derivatives of inputs.
- `mtkcompile=nothing`: Controls the simplification backend. `false` (default) uses NetworkDynamics'
  own algebraic reduction pipeline. `true` delegates to MTK's `mtkcompile`: when `ModelingToolkit` is
  loaded this performs full structural simplification; when only `ModelingToolkitBase` is available it
  falls back to the stub `mtkcompile` (no structural simplification). `nothing` inherits the global
  default set via `NetworkDynamics.set_mtkcompile!`. `:compare` runs both and prints a comparison
  (for debugging).
"""
function EdgeModel(
    sys::System, srcin, dstin, srcout, dstout;
    verbose=false,
    name=getname(sys),
    ff_to_constraint=false,
    assume_io_coupling=false,
    extin=nothing,
    mtkcompile=nothing,
    kwargs...
)
    mtkcompile = isnothing(mtkcompile) ? MTKCOMPILE_DEFAULT[] : mtkcompile
    if mtkcompile==:compare
        args = (sys, srcin, dstin, srcout, dstout)
        kwargs = (; verbose, name, ff_to_constraint, assume_io_coupling, extin, kwargs...)
        return _compare_mtkcompile(EdgeModel, args, kwargs)
    end

    warn_missing_features(sys)
    sys = resolve_bound_to(sys)
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

    gen = generate_io_function_cached(sys, ins, outs; verbose, ff_to_constraint, assume_io_coupling, mtkcompile)

    f = gen.f
    g = singlesided ? gwrap(gen.g) : gen.g
    obsf = gen.obsf

    sysdefaults = initial_conditions(sys)
    sysguesses = ModelingToolkitBase.guesses(sys)
    _sym = getname.(gen.states)
    sym = [s => _get_metadata(sys, s, sysdefaults, sysguesses) for s in _sym]

    _psym = getname.(gen.params)
    psym = map(gen.params) do p
        pname = getname(p)
        md = _get_metadata(sys, pname, sysdefaults, sysguesses)
        if p in gen.unused_params
            md[:unused] = true
        end
        pname => md
    end

    _obssym = getname.(gen.obsstates)
    obssym = [s => _get_metadata(sys, s, sysdefaults, sysguesses) for s in _obssym]

    _insym_src = getname.(srcin)
    insym_src = [s => _get_metadata(sys, s, sysdefaults, sysguesses)  for s in _insym_src]
    _insym_dst = getname.(dstin)
    insym_dst = [s => _get_metadata(sys, s, sysdefaults, sysguesses)  for s in _insym_dst]
    insym = (;src=insym_src, dst=insym_dst)

    if singlesided
        _outsym_dst = getname.(dstout)
        outsym = [s => _get_metadata(sys, s, sysdefaults, sysguesses)  for s in _outsym_dst]
    else
        _outsym_src = getname.(srcout)
        outsym_src = [s => _get_metadata(sys, s, sysdefaults, sysguesses)  for s in _outsym_src]
        _outsym_dst = getname.(dstout)
        outsym_dst = [s => _get_metadata(sys, s, sysdefaults, sysguesses)  for s in _outsym_dst]
        outsym = (;src=outsym_src, dst=outsym_dst)
    end

    mass_matrix = gen.mass_matrix
    c = EdgeModel(;f, g, sym, insym, outsym, psym, obssym,
            obsf, mass_matrix, ff=gen.fftype, name, extin=extin_nwidx,
            aliasmap=gen.aliasmap, allow_output_sym_clash=true, kwargs...)
    set_metadata!(c, :observed, gen.observed)
    set_obsrules!(c, gen.obsrules)
    set_metadata!(c, :equations, gen.equations)
    set_metadata!(c, :full_equations, gen.full_equations)
    set_metadata!(c, :outputeqs, gen.outputeqs)
    set_metadata!(c, :full_outputeqs, gen.full_outputeqs)
    set_metadata!(c, :odesystem, sys)

    for formula in something(gen.initformulas, ())
        add_initformula_lenient!(c, formula)
    end
    for formula in something(gen.guessformulas, ())
        add_guessformula_lenient!(c, formula)
    end

    # apply postprocessing functions from subcomponent metadata
    apply_component_postprocessing!(c)
    c
end

"""
For a given system and name, extract all the relevant meta we want to keep for the component model.

since defaults(sys) and guesses(sys) is relativly expensive those need to be provided
"""
function _get_metadata(sys, name, alldefaults, allguesses)
    md = Dict{Symbol,Any}()
    sym = try
        getproperty_symbolic(sys, name; might_contain_toplevel_ns=false)
    catch e
        if !endswith(string(name), "ˍt") # known for "internal" derivatives
            @warn "Could not extract metadata for $name $(e.msg)"
        end
        return md
    end
    if haskey(alldefaults, sym)
        def = unwrap_const(alldefaults[sym])
        md[:default] = def
    end

    if haskey(allguesses, sym)
        guess = unwrap_const(allguesses[sym])
        if is_symbolic_constant(guess)
            md[:guess] = symbolic_constant_value(guess)
        else
            throw(ArgumentError(
                "Symbolic guess `$name => $guess` is not supported: the `guess` option only \
                 takes constant values. Use the `guessf` variable option \
                 (`@variables $name [guessf = $guess]`) or `set_guessf` to declare a guess \
                 formula instead."))
        end
    end

    if ModelingToolkitBase.hasbounds(sym)
        md[:bounds] = ModelingToolkitBase.getbounds(sym)
    end
    if ModelingToolkitBase.hasdescription(sym)
        md[:description] = ModelingToolkitBase.getdescription(sym)
    end

    default_from = Symbolics.getmetadata(unwrap(sym), ParameterDefaultFrom, nothing)
    if !isnothing(default_from)
        # arrives as a `(:src, :busbar₊S_b)` tuple of QuoteNodes from @mtkmodel/@component
        md[:default_from] = normalize_default_from(default_from)
    end
    md
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
                              expression=Val{false}, verbose=false, mtkcompile=false,
                              ff_to_constraint, assume_io_coupling)
    # TODO: scalarize vector symbolics/equations?

    # f_* may be given in namepsace version or as symbols, resolve to unnamespaced Symbolic
    inputss = map(inputss) do in
        getproperty_symbolic.(Ref(_sys), in)
    end
    allinputs = convert(Vector{ST}, reduce(union, inputss))
    outputss = map(outputss) do out
        getproperty_symbolic.(Ref(_sys), out)
    end
    alloutputs = convert(Vector{ST}, reduce(union, outputss))

    # Both must see the hierarchical system: `initf` metadata expressions are only
    # namespaceable before flattening, and mtkcompile synthesizes bindings of its own.
    assert_no_state_bindings(_sys)
    initf = collect_initf(_sys)
    guessf = collect_guessf(_sys)

    # always expand connections before simplification
    _sys = ModelingToolkitBase.expand_connections(_sys)
    _ivcand = independent_variables(_sys)
    if length(_ivcand) == 1
        iv = only(independent_variables(_sys))
        # lots of places just assume time is called `:t` instead of threading the name around,
        # so assert it once here
        getname(iv) === :t || throw(ArgumentError(
            "The independent variable of $(nameof(_sys)) is :$(getname(iv)); NetworkDynamics \
            requires it to be named `t`."))
    elseif length(_ivcand) > 1
        throw(ArgumentError(
            "The system $(nameof(_sys)) has multiple independent variables: \
            $(getname.(_ivcand)). NetworkDynamics only supports a single independent variable."))
    else
        # its probably t
        iv = ModelingToolkitBase.t_nounits
    end

    # assume_io_coupling means, we expect a direct dependency chain output -> input
    # we fake this, by replacing all inputs with `input + implicit_output(outputs...)`
    # which forces MTK to consider the dependency
    if assume_io_coupling
        impl_outputs = implicit_output(sum(alloutputs))
        subs = Dict(in => in + impl_outputs for in in allinputs)
        _expanded_eqs = ModelingToolkitBase.get_eqs(_sys)
        for i in eachindex(_expanded_eqs)
            _expanded_eqs[i] = substitute(_expanded_eqs[i], subs)
        end
    end

    eqs, obseqs_sorted, states, params = if mtkcompile
        simplify_with_mtkcompile(_sys, allinputs, alloutputs; verbose)
    else
        simplify_without_mtkcompile(_sys, allinputs, alloutputs; verbose, ff_to_constraint)
    end

    # get rid of the implicit_output(⋅) terms
    remove_implicit_output_fn!(eqs)
    remove_implicit_output_fn!(obseqs_sorted)

    eqs, obseqs_sorted, states, aliasmap = pick_best_alias_names(eqs, obseqs_sorted, states, alloutputs, allinputs; verbose)

    # check that there are no rhs differentials in the equations
    if !isempty(rhs_differentials(vcat(eqs, obseqs_sorted)))
        diffs = rhs_differentials(vcat(eqs, obseqs_sorted))
        throw(RHSDifferentialsError(repr.(diffs)))
    end

    # obs can only depend on parameters (including allinputs) or states
    obs_subs = OrderedDict(eq.lhs => eq.rhs for eq in obseqs_sorted)
    obs_deps = setdiff(_all_rhs_symbols(obs_subs), Set(keys(obs_subs))) # do not count self-dependency
    if !(obs_deps ⊆ Set(params) ∪ Set(allinputs) ∪ Set(states) ∪ iv)
        @warn "obs_deps !⊆ params ∪ inputs ∪ unknowns. Difference: $(setdiff(obs_deps, Set(params) ∪ Set(allinputs) ∪ Set(states)))"
    end

    # find the output equations, this might remove them from obseqs_sorted (obs_subs stays intact)
    outeqs = Equation[]
    for out in Iterators.flatten(outputss)
        if out ∈ Set(states)
            push!(outeqs, out ~ out)
        elseif out ∈ keys(obs_subs)
            # if its a observed, we need to check for ff behavior
            fulldeps = _all_dependencies(obs_subs[out], obs_subs);
            if ff_to_constraint && !isempty(fulldeps ∩ allinputs)
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
        for (i, eq) in enumerate(eqs)
           if eq_type(eq)[1] == :explicit_algebraic
               error("Can't have explicit algebraic equations at this point")
           end
        end
        # create massmatrix, we don't use the method provided by System because of reordering
        mm = generate_massmatrix(eqs)
        verbose && @info "Generated mass matrix" mm
        mm
    end

    out_deps = _all_dependencies(outeqs, obs_subs)
    fftype = _determine_fftype(out_deps, states, allinputs, params, iv)

    # filter out unnecessary parameters
    unused_params = let
        # we need to collect obs and var deps separately, because replacement of obs in vars might lead to
        # symbolic simplifications which we don't have in the actual equations later!
        # i.e. dt(x) ~ x - x0 and x ~ x0 + y leads to dt(x) ~ y in substitution but is not resolved in actual f
        var_deps = _all_rhs_symbols(eqs)
        Set(setdiff(params, (var_deps ∪ obs_deps ∪ out_deps))) # do not exclude obs_deps
    end
    if verbose && !isempty(unused_params)
        @info "Parameters $(unused_params) do not appear in equations of f and g and will be marked as unused."
    end

    # now generate the actual functions
    if !isempty(eqs)
        f_ip = build_fg(eqs, obs_subs, (states, inputss..., params, iv); expression)
    else
        f_ip = nothing
    end

    gformargs = if fftype isa PureFeedForward
        (inputss..., params, iv)
    elseif fftype isa FeedForward
        (states, inputss..., params, iv)
    elseif fftype isa NoFeedForward
        (states, params, iv)
    elseif fftype isa PureStateMap
        (states,)
    end
    _g_ip = build_fg(outeqs, obs_subs, gformargs; expression)
    # for more than 1 output, wrap function in MultipleOutputWrapper
    g_ip = if length(outputss) == 1
        _g_ip
    else
        MultipleOutputWrapper{fftype, length(outputss), typeof(_g_ip)}(_g_ip)
    end

    obsstates = [eq.lhs for eq in obseqs_sorted]
    if !isempty(obsstates)
        obsf_ip = build_obsf(obsstates, obs_subs, allinputs, (states, inputss..., params, iv); expression)
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
        full_equations = [eq.lhs ~ fixpoint_sub(eq.rhs, obs_subs) for eq in eqs],
        full_outputeqs = [eq.lhs ~ fixpoint_sub(eq.rhs, obs_subs) for eq in outeqs],
        observed=obseqs_sorted,
        obsrules = build_obsrules(obseqs_sorted, outeqs, eqs, aliasmap, iv),
        aliasmap,
        iv,
        params,
        unused_params,
        initformulas = initf_to_initformulas(initf),
        guessformulas = guessf_to_guessformulas(guessf)
    )
end

function _collect_deps_on_obs(terms, obs_subs)
    deps = Set{ST}()
    for term in terms
        _collect_deps_on_obs!(deps, obs_subs, term)
    end
    deps
end
function _collect_deps_on_obs!(deps, obs_subs, term)
    termdeps = get_variables_deriv(term)
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

_all_rhs_symbols(term) = get_variables_deriv(term)
_all_rhs_symbols(eq::Equation) = get_variables_deriv(eq.rhs)
_all_rhs_symbols(eqs::Union{AbstractVector,AbstractDict}) = mapreduce(eq->get_variables_deriv(eq isa Pair ? eq.second : eq.rhs), ∪, eqs, init=Set{ST}())

"""
Search for recursive dependencies in `term` given a dictionary `dict` of substitutions.
"""
function _all_dependencies(term, dict)
    deps = Set{ST}()
    _recursive_collect_dependencies!(deps, term, dict)
    deps
end
function _recursive_collect_dependencies!(deps, term, dict)
    termdeps = _all_rhs_symbols(term)
    for sym in termdeps
        if haskey(dict, sym)
            _recursive_collect_dependencies!(deps, dict[sym], dict)
        else
            push!(deps, sym)
        end
    end
    deps
end

using PrecompileTools: @compile_workload
@compile_workload begin
    # include("MTKExt_precomp_workload.jl")
end

end
