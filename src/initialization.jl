"""
    find_fixpoint(nw::Network, [x0::NWState=NWState(nw)], [p::NWParameter=x0.p]; kwargs...)
    find_fixpoint(nw::Network, x0::AbstractVector, p::AbstractVector; kwargs...)
    find_fixpoint(nw::Network, x0::AbstractVector; kwargs...)

Convenience wrapper around `SteadyStateProblem` from SciML-ecosystem.
Constructs and solves the steady state problem, returns found value wrapped as `NWState`.
"""
function find_fixpoint(nw::Network, x0::AbstractVector; kwargs...)
    find_fixpoint(nw, x0, NWParameter(nw); kwargs...)
end
function find_fixpoint(nw::Network, p::NWParameter; kwargs...)
    find_fixpoint(nw, NWState(p; ufill=0), p; kwargs...)
end
function find_fixpoint(nw::Network,
                       x0::NWState=NWState(nw; ufill=0),
                       p::NWParameter=x0.p;
                       kwargs...)
    find_fixpoint(nw, uflat(x0), pflat(p); kwargs...)
end
function find_fixpoint(nw::Network, x0::AbstractVector, p::AbstractVector;
                       alg=SSRootfind(), kwargs...)
    if any(isnan, x0) || any(isnan, p)
        missing_vars = String[]
        missing_params = String[]

        if any(isnan, x0)
            nan_indices = findall(isnan, x0)
            variable_symbols = NetworkDynamics.SII.variable_symbols(nw)
            missing_vars = [string(variable_symbols[i]) for i in nan_indices]
        end

        if any(isnan, p)
            nan_indices = findall(isnan, p)
            parameter_symbols = NetworkDynamics.SII.parameter_symbols(nw)
            missing_params = [string(parameter_symbols[i]) for i in nan_indices]
        end

        if !isempty(missing_vars) || !isempty(missing_params)
            error_msg = "find_fixpoint inputs contain NaNs, indicating missing default values:\n"
            if !isempty(missing_vars)
                error_msg *= "  Variables:  " * join(missing_vars, ", ") * "\n"
            end
            if !isempty(missing_params)
                error_msg *= "  Parameters: " * join(missing_params, ", ") * "\n"
            end
            error_msg *= "This may indicate missing default values in component metadata or explicitly passed inputs with NaNs."
            throw(ArgumentError(error_msg))
        end
    end
    prob = SteadyStateProblem(nw, x0, p)
    sol = _solve_fixpoint(prob, alg; kwargs...)
    if !SciMLBase.successful_retcode(sol.retcode)
        @warn "Solver did not finish retcode = $(sol.retcode) (alg = $alg)!"
        error("""
        Could not find fixpoint, solver returned $(sol.retcode)! For debugging, \
        it is advised to manually construct the steady state problem and try \
        different solvers/arguments:

        prob = SteadyStateproblem(nw, uflat(nwstat), pflat(nwpara))
        sol = solve(prob, alg; kwargs...)
        x0 = NWState(nw, sol.u)

        For detail see https://docs.sciml.ai/NonlinearSolve/stable/native/steadystatediffeq/
        """)
    end
    NWState(nw, sol.u, NWParameter(nw, p))
end

function _solve_fixpoint(prob, alg::AbstractNonlinearSolveAlgorithm; kwargs...)
    _prob = NonlinearProblem(prob)
    SciMLBase.solve(_prob, alg; kwargs...)
end

function _solve_fixpoint(prob, alg::SteadyStateDiffEqAlgorithm; kwargs...)
    sol = SciMLBase.solve(prob, alg; kwargs...)
end

function initialization_problem(cf::T,
    defaults,
    guesses,
    bounds,
    initconstraint=nothing;
    t=NaN,
    apply_bound_transformation=true,
    verbose=true
) where {T<:ComponentModel}
    hasinsym(cf) || throw(ArgumentError("Component model musst have `insym`!"))

    # The _m[s] suffix means a bitmask which indicate the free variables
    # the _fix[s] suffix means an array of same length as the symbols, which contains the fixed values
    # the [s] only appears for multiple masks, i.e. (out1, out2, ..) and (in1, in2, ...)

    outfree_ms = Tuple(
        [!haskey(defaults, s) for s in sv]
        for sv in outsym_normalized(cf)
    )
    outfixs = Tuple(
        Float64[haskey(defaults, s) ? defaults[s] : NaN for s in sv]
        for sv in outsym_normalized(cf)
    )

    ufree_m = [!haskey(defaults, s) for s in sym(cf)]
    ufix = Float64[haskey(defaults, s) ? defaults[s] : NaN for s in sym(cf)]

    infree_ms = Tuple(
        [!haskey(defaults, s) for s in sv]
        for sv in insym_normalized(cf)
    )
    infixs = Tuple(
        Float64[haskey(defaults, s) ? defaults[s] : NaN for s in sv]
        for sv in insym_normalized(cf)
    )

    is_free_p = s -> !is_unused(cf, s) && !haskey(defaults, s) # free p are not unused and have no default
    pfree_m = [is_free_p(s) for s in psym(cf)]
    pfix = Float64[haskey(defaults, s) ? defaults[s] : NaN for s in psym(cf)]

    # count free variables and equations
    Nfree = mapreduce(sum, +, outfree_ms) + sum(ufree_m) + mapreduce(sum, +, infree_ms) + sum(pfree_m)

    # return "fake NonlinearProblem" if no free variables
    iszero(Nfree) && return ((; u0=[]), identity)

    # check for additonal equations
    additional_Neqs = isnothing(initconstraint) ? 0 : dim(initconstraint)

    Neqs = dim(cf) + mapreduce(length, +, outfree_ms) + additional_Neqs

    freesym = vcat(mapreduce((syms, map) -> syms[map], vcat, outsym_normalized(cf), outfree_ms),
                   sym(cf)[ufree_m],
                   mapreduce((syms, map) -> syms[map], vcat, insym_normalized(cf), infree_ms),
                   psym(cf)[pfree_m])

    @assert length(freesym) == Nfree
    if Neqs < Nfree
        throw(
            ArgumentError("Initialization problem underconstraint. \
                           $(Neqs) Equations for $(Nfree) free variables: $freesym. Consider \
                           passing additional constraints using `InitConstraint`.")
        )

    end

    # which parts of the nonlinear u (unl) correspond to which free parameters
    nextrange = let lasti = 0
        (mask) -> begin
            start = lasti + 1
            stop = lasti + sum(mask)
            lasti = stop
            start:stop
        end
    end
    unl_range_outs = map(nextrange, outfree_ms)
    unl_range_u = nextrange(ufree_m)
    unl_range_ins = map(nextrange, infree_ms)
    unl_range_p = nextrange(pfree_m)
    @assert vcat(unl_range_outs..., unl_range_u, unl_range_ins..., unl_range_p) == 1:Nfree

    # check for positivity and negativity constraints
    bound_types = map(freesym) do sym
        if haskey(bounds, sym)
            bound = bounds[sym]
            if bound[1] >= 0 && bound[2] > bound[1]
                return :pos
            elseif bound[1] < bound[2]  && bound[2] <= 0
                return :neg
            end
        end
        :none
    end
    if !apply_bound_transformation || all(isequal(:none), bound_types)
        boundT! = identity
        inv_boundT! = identity
    else
        if verbose
            idxs = findall(!isequal(:none), bound_types)
            @info "Apply positivity/negativity conserving variable transformation on $(freesym[idxs]) to satisfy bounds."
        end
        boundT! = (u) -> begin
            for i in eachindex(u, bound_types)
                if bound_types[i] == :pos
                    u[i] = u[i]^2
                elseif bound_types[i] == :neg
                    u[i] = -u[i]^2
                end
            end
            return u
        end
        inv_boundT! = (u) -> begin
            for i in eachindex(u, bound_types)
                if bound_types[i] == :pos
                    u[i] = sqrt(u[i])
                elseif bound_types[i] == :neg
                    u[i] = sqrt(-u[i])
                end
            end
            return u
        end
    end

    # check for missing guesses
    missing_guesses = Symbol[]
    uguess = map(freesym) do s
        if haskey(guesses, s)
            Float64(guesses[s])
        else
            push!(missing_guesses, s)
        end
    end
    isempty(missing_guesses) || throw(ArgumentError("Missing guesses for free variables $(missing_guesses)"))

    # apply bound conserving transformation to initial state
    try
        inv_boundT!(uguess)
    catch e
        throw(ArgumentError("Failed to apply bound-conserving transformation to initial guess. \
                             This typically happens when a guess has the wrong sign for its bounds \
                             (e.g., negative guess for a positively-bounded variable). Original error: $(e)"))
    end

    chunksize = ForwardDiff.pickchunksize(Nfree)
    fg = compfg(cf)
    unlcache = map(d->DiffCache(zeros(d), chunksize), length(freesym))
    outcaches = map(d->DiffCache(zeros(d), chunksize), outdim_normalized(cf))
    ucache = DiffCache(zeros(dim(cf)), chunksize)
    incaches = map(d->DiffCache(zeros(d), chunksize), indim_normalized(cf))
    pcache = DiffCache(zeros(pdim(cf)), chunksize)

    # generate a function to apply the additional constraints
    additional_cf = prep_initiconstraint(cf, initconstraint, chunksize) #is noop for add_c == nothing

    fz = (dunl, unl, _) -> begin
        # apply the bound conserving transformation
        unlbuf = PreallocationTools.get_tmp(unlcache, unl)
        unlbuf .= unl
        boundT!(unlbuf)

        outbufs = PreallocationTools.get_tmp.(outcaches, Ref(dunl))
        ubuf = PreallocationTools.get_tmp(ucache, dunl)
        inbufs = PreallocationTools.get_tmp.(incaches, Ref(dunl))
        pbuf = PreallocationTools.get_tmp(pcache, dunl)

        # prefill buffers with fixed values
        for (buf, fix) in zip(outbufs, outfixs)
            buf .= fix
        end
        ubuf .= ufix
        for (buf, fix) in zip(inbufs, infixs)
            buf .= fix
        end
        pbuf .= pfix

        # overwrite nonfixed values
        for (buf, mask, range) in zip(outbufs, outfree_ms, unl_range_outs)
            _overwrite_at_mask!(buf, mask, unlbuf, range)
        end
        _overwrite_at_mask!(ubuf, ufree_m, unlbuf, unl_range_u)
        for (buf, mask, range) in zip(inbufs, infree_ms, unl_range_ins)
            _overwrite_at_mask!(buf, mask, unlbuf, range)
        end
        _overwrite_at_mask!(pbuf, pfree_m, unlbuf, unl_range_p)

        # view into du buffer for the fg funtion
        @views dunl_fg = dunl[1:dim(cf)]
        # view into the output buffer for the outputs
        @views dunl_out = dunl[dim(cf)+1:dim(cf)+reduce(+,outdim(cf))]
        # view into the output buffer to the additional constraints
        @views dunl_additional = dunl[dim(cf)+reduce(+,outdim(cf))+1:end]

        # this fills the second half of the du buffer with the fixed and current outputs
        dunl_out .= RecursiveArrayTools.ArrayPartition(outbufs...)
        # execute fg to fill dunl and outputs
        fg(outbufs, dunl_fg, ubuf, inbufs, pbuf, t)
        # calculate the residual for the second half ov the dunl buf, the outputs
        dunl_out .= dunl_out .- RecursiveArrayTools.ArrayPartition(outbufs...)
        # execute the additonal constraints
        additional_cf(dunl_additional, outbufs, ubuf, inbufs, pbuf, t)
        nothing
    end

    nlf = NonlinearFunction(fz; resid_prototype=zeros(Neqs), sys=SII.SymbolCache(freesym))

    if Neqs == Nfree
        verbose && @info "Initialization problem is fully constrained. Created NonlinearLeastSquaresProblem for $freesym"
    else
        verbose && @info "Initialization problem is overconstrained ($Nfree vars for $Neqs equations). Create NonlinearLeastSquaresProblem for $freesym."
    end
    (NonlinearLeastSquaresProblem(nlf, uguess), boundT!)
end
function _overwrite_at_mask!(target, mask, source, range)
    src_v = view(source, range)
    j = 1
    for i in eachindex(target, mask)
        if mask[i]
            target[i] = src_v[j]
            j += 1
        end
    end
end

"""
    initialize_component(cf;
                         defaults=get_defaults_dict(cf),
                         guesses=get_guesses_dict(cf),
                         bounds=get_bounds_dict(cf),
                         default_overrides=nothing,
                         guess_overrides=nothing,
                         bound_overrides=nothing,
                         additional_initformula=nothing,
                         additional_initconstraint=nothing,
                         verbose=true,
                         apply_bound_transformation=true,
                         t=NaN,
                         tol=1e-10,
                         kwargs...)

The function solves a nonlinear problem to find values for all free variables/parameters
(those without defaults) that satisfy the component equations in steady state
(i.e. RHS equals 0). The initial guess for each variable depends on the provided
`guesses` parameter (defaults to the metadata `guess` values).

## Parameters
- `cf`: ComponentModel to initialize
- `defaults`: Dictionary of default values (defaults to metadata defaults)
- `guesses`: Dictionary of initial guesses (defaults to metadata guesses)
- `bounds`: Dictionary of bounds (defaults to metadata bounds)
- `default/guess/bound_overrides`: Dictionary to merge with `defaults`/`guesses`/`bounds`.
  You can use `nothing` as a value for any key to remove that entry from the respective dictionary.
- `additional_initformula`: Additional initialization formulas to apply beyond those in component metadata
- `additional_initconstraint`: Additional initialization constraints to apply beyond those in component metadata
- `verbose`: Whether to print information during initialization
- `apply_bound_transformation`: Whether to apply bound-conserving transformations
- `t`: Time at which to solve for steady state. Only relevant for components with explicit time dependency.
- `tol`: Tolerance for the residual of the initialized model (defaults to `1e-10`). Init throws error if resid < tol.
- `kwargs...`: Additional arguments passed to the nonlinear solver

## Returns
- Dictionary mapping symbols to their values (complete state including defaults and initialized values)

## Bounds of free variables
When encountering any bounds in the free variables, NetworkDynamics will try to conserve them
by applying a coordinate transformation. This behavior can be suppressed by setting `apply_bound_transformation=false`.
The following transformations are used:
- (a, b) intervals where both a and b are positive are transformed to `u^2`/`sqrt(u)`
- (a, b) intervals where both a and b are negative are transformed to `-u^2`/`sqrt(-u)`
"""
function initialize_component(cf;
                             defaults=get_defaults_dict(cf),
                             guesses=get_guesses_dict(cf),
                             bounds=get_bounds_dict(cf),
                             default_overrides=nothing,
                             guess_overrides=nothing,
                             bound_overrides=nothing,
                             additional_initformula=nothing,
                             additional_initconstraint=nothing,
                             verbose=true,
                             apply_bound_transformation=true,
                             t=NaN,
                             tol=1e-10,
                             kwargs...)

    defaults = isnothing(default_overrides) ? defaults : merge(defaults, default_overrides)
    guesses  = isnothing(guess_overrides)   ? guesses  : merge(guesses, guess_overrides)
    bounds   = isnothing(bound_overrides)   ? bounds   : merge(bounds, bound_overrides)

    # filter out nothing values
    defaults = filter(p -> !isnothing(p.second), defaults)
    guesses = filter(p -> !isnothing(p.second), guesses)
    bounds = filter(p -> !isnothing(p.second), bounds)

    # Extract metadata and merge with additional constraints/formulas
    metadata_initformulas = has_initformula(cf) ? get_initformulas(cf) : nothing
    combined_initformulas = collect_initformulas(metadata_initformulas, additional_initformula)

    # Apply initialization formulas to defaults
    if !isnothing(combined_initformulas)
        apply_init_formulas!(defaults, combined_initformulas; verbose)
    end

    metadata_constraint = has_initconstraint(cf) ? get_initconstraints(cf) : nothing
    combined_constraint = merge_initconstraints(metadata_constraint, additional_initconstraint)

    prob, boundT! = initialization_problem(cf, defaults, guesses, bounds,
                                          combined_constraint;
                                          verbose, apply_bound_transformation, t)

    # Dictionary for complete state (defaults + initialized values)
    init_state = Dict{Symbol, Float64}()
    observable_defaults = Dict{Symbol, Float64}()
    _obssym = obssym(cf)
    for (sym, val) in defaults
        if sym in _obssym
            observable_defaults[sym] = val
        else
            init_state[sym] = val
        end
    end

    if !isempty(prob.u0)
        sol = SciMLBase.solve(prob; verbose, kwargs...)

        if sol.prob isa NonlinearLeastSquaresProblem && sol.retcode == SciMLBase.ReturnCode.Stalled
            res = LinearAlgebra.norm(sol.resid)
            @warn "Initialization for component stalled with residual $(res)"
        elseif !SciMLBase.successful_retcode(sol.retcode)
            throw(ArgumentError("Initialization failed. Solver returned $(sol.retcode)"))
        else
            res = LinearAlgebra.norm(sol.resid)
            verbose && @info "Initialization successful with residual $(res)"
        end

        # Transform back to original space
        u = boundT!(copy(sol.u))

        # Add solved values to the complete dictionary
        free_symbols = SII.variable_symbols(sol)
        for (sym, val) in zip(free_symbols, u)
            init_state[sym] = val
        end
    else
        res = init_residual(cf, init_state, t=NaN)
        verbose && @info "No free variables! Residual $(res)"
    end
    if !(res < tol)
        error("Initialized model has a residual larger then specified tolerance $(res) > $(tol)! \
               Fix initialization or increase tolerance to supress error.")
    end


    # Check for broken bounds using the complete state
    broken_bnds = broken_bounds(cf, init_state, bounds)
    if !isempty(broken_bnds)
        broken_msgs = ["$sym = $val (bounds: $lb..$ub)" for (sym, val, (lb, ub)) in broken_bnds]
        @warn "Initialized model has broken bounds. Try to adapt the initial guesses!" *
              "\n" * join(broken_msgs, "\n")
    end

    # Check for broken observable defaults
    if !isempty(observable_defaults)
        broken_obs = broken_observable_defaults(cf, init_state, observable_defaults)
        if !isempty(broken_obs)
            broken_msgs = ["$sym = $val (default: $def)" for (sym, def, val) in broken_obs]
            @warn "Initialized model has observables that differ from their specified defaults:" *
                  "\n" * join(broken_msgs, "\n")
        end
    end

    return init_state
end

"""
    initialize_component!(cf::ComponentModel;
                          defaults=nothing,
                          guesses=nothing,
                          bounds=nothing,
                          default_overrides=nothing,
                          guess_overrides=nothing,
                          bound_overrides=nothing,
                          additional_initformula=nothing,
                          additional_initconstraint=nothing,
                          verbose=true,
                          t=NaN,
                          kwargs...)

Mutating version of [`initialize_component`](@ref). See this docstring for all details.
In contrast to the non mutating version, this function reads in defaults and guesses from the
symbolic metadata and writes the initialized values back in to the metadata.

## Parameters
- `cf`: ComponentModel to initialize
- `defaults`: Optional dictionary to replace all metadata defaults
- `guesses`: Optional dictionary to replace all metadata guesses
- `bounds`: Optional dictionary to replace all metadata bounds
- `default/guess/bound_overrides`: Dict of values that override existing
   default/guess/bound metadata. Use `nothing` as a value for any key to remove
   that metadata entry from the component model.
- `additional_initformula`: Additional initialization formulas to apply beyond those in component metadata
- `additional_initconstraint`: Additional initialization constraints to apply beyond those in component metadata
- `verbose`: Whether to print information during initialization
- `t`: Time at which to solve for steady state. Only relevant for components with explicit time dependency.
- All other `kwargs` are passed to `initialize_component`

When `defaults`, `guesses`, or `bounds` are provided, they replace the corresponding
metadata in the component model. Any keys in the original metadata that are not in
the provided dictionaries will be removed, and new keys will be added.
"""
function initialize_component!(cf;
                              defaults=nothing,
                              guesses=nothing,
                              bounds=nothing,
                              default_overrides=nothing,
                              guess_overrides=nothing,
                              bound_overrides=nothing,
                              additional_initformula=nothing,
                              additional_initconstraint=nothing,
                              verbose=true,
                              t=NaN,
                              kwargs...)

    # Synchronize defaults, guesses, and bounds
    _sync_metadata!(cf, defaults, get_defaults_dict, set_default!, delete_default!, "default", verbose)
    _sync_metadata!(cf, guesses,  get_guesses_dict,  set_guess!,   delete_guess!,   "guess",   verbose)
    _sync_metadata!(cf, bounds,   get_bounds_dict,   set_bounds!,  delete_bounds!,  "bounds",  verbose)

    for (name, set_fn!, rm_fn!, dict) in (
        ("default", set_default!, delete_default!, default_overrides),
        ("guess", set_guess!, delete_guess!, guess_overrides),
        ("bound", set_bounds!, delete_bounds!, bound_overrides)
    )
        isnothing(dict) && continue
        for (sym, val) in dict
            if !isnothing(val)
                set_fn!(cf, sym, val)
                verbose && println("Set additional $name for $sym: $val")
            else
                rm_fn!(cf, sym)
                verbose && println("Remove $name for $sym")
            end
        end
    end

    # Now proceed with initialization using the updated metadata
    init_state = initialize_component(
        cf;
        additional_initformula=additional_initformula,
        additional_initconstraint=additional_initconstraint,
        verbose=verbose,
        t=t,
        kwargs...  # Only pass the remaining kwargs
    )

    # Calculate residual for validation
    resid = init_residual(cf, init_state; t=t)

    # delete all inits
    for s in keys(get_inits_dict(cf))
        delete_init!(cf, s)
    end

    # Write back the initialized values to the component metadata
    for (sym, val) in init_state
        if has_default(cf, sym)
            # Check if initialized value differs from existing default
            # This can happen when init formulas override default values
            get_default(cf, sym) != val && set_default!(cf, sym, val)
        else
            # No default exists, store as init value
            set_init!(cf, sym, val)
        end
    end

    cf
end
function _sync_metadata!(cf, provided, get_orig_fn, set_fn!, remove_fn!, type_name, verbose)
    isnothing(provided) && return

    original = get_orig_fn(cf)

    # Identify keys to be added/updated/removed
    provided_keys = Set(keys(provided))
    original_keys = Set(keys(original))

    missing_keys = setdiff(original_keys, provided_keys)
    new_keys = setdiff(provided_keys, original_keys)
    common_keys = intersect(provided_keys, original_keys)
    changed_keys = filter(k -> provided[k] != original[k], collect(common_keys))

    # Remove missing keys
    for sym in missing_keys
        remove_fn!(cf, sym)
        verbose && println("Removed $type_name for $sym (was $(original[sym]))")
    end

    # Add new keys
    for sym in new_keys
        set_fn!(cf, sym, provided[sym])
        verbose && println("Added $type_name for $sym: $(provided[sym])")
    end

    # Update changed keys
    for sym in changed_keys
        set_fn!(cf, sym, provided[sym])
        verbose && println("Updated $type_name for $sym: $(original[sym]) → $(provided[sym])")
    end
end


function isinitialized(cf::ComponentModel)
    all(has_default_or_init(cf, s) || is_unused(cf, s) for s in vcat(sym(cf), psym(cf)))
end

"""
    init_residual(cf::ComponentModel, [state=get_defaults_or_inits_dict(cf)]; t=NaN)

Calculates the residual |du| for the given component model using the values provided.
If no state dictionary is provided, it uses the values from default/init [Metadata](@ref).

See also [`initialize_component`](@ref).
"""
function init_residual(cf::ComponentModel, state=get_defaults_or_inits_dict(cf); t=NaN)
    # Check that all necessary symbols are present in the state dictionary
    needed_symbols = Set(vcat(
        sym(cf),                                # states
        filter(s->!is_unused(cf, s), psym(cf)), # used parameters
        insym_flat(cf),                         # inputs
        outsym_flat(cf)                         # outputs
    ))

    missing = setdiff(needed_symbols, keys(state))
    if !isempty(missing)
        throw(ArgumentError("State dictionary is missing required symbols: $missing. \
                             These symbols need values to calculate the residual."))
    end

    outs = Tuple(Float64[get(state, s, NaN) for s in sv] for sv in outsym_normalized(cf))
    u = Float64[get(state, s, NaN) for s in sym(cf)]
    ins = Tuple(Float64[get(state, s, NaN) for s in sv] for sv in insym_normalized(cf))
    p = Float64[is_unused(cf, s) ? NaN : get(state, s, NaN) for s in psym(cf)]

    # collect additional constraints
    additional_constraint = if has_initconstraint(cf)
        _c = get_initconstraints(cf)
        length(_c) == 1 ? only(_c) : InitConstraint(_c...)
    else
        nothing
    end
    additional_Neqs = isnothing(additional_constraint) ? 0 : dim(additional_constraint)
    additional_cf = prep_initiconstraint(cf, additional_constraint, 0) #is noop for add_c == nothing

    Nout = reduce(+, outdim(cf))
    res = zeros(dim(cf) + Nout + additional_Neqs)

    res_fg = @views res[1:dim(cf)]
    res_out = @views res[dim(cf)+1:dim(cf)+Nout]
    res_add = @views res[dim(cf)+Nout+1:end]

    res_out .= RecursiveArrayTools.ArrayPartition(outs...) # fill with provided outputs
    compfg(cf)(outs, res_fg, u, ins, p, t)
    res_out .= res_out .- RecursiveArrayTools.ArrayPartition(outs...) # compare to calculated ouputs

    additional_cf(res_add, outs, u, ins, p, t) # apply additional constraints

    return LinearAlgebra.norm(res)
end

function broken_bounds(cf, state=get_defaults_or_inits_dict(cf), bounds=get_bounds_dict(cf))
    vals = get_initial_state(cf, state, keys(bounds))

    broken = []
    for (val, sym, bound) in zip(vals, keys(bounds), values(bounds))
        if !bounds_satisfied(val, bound)
            push!(broken, (sym, val, bound))
        end
    end
    broken
end
function bounds_satisfied(val, bounds)
    @assert length(bounds) == 2
    !isnothing(val) && !isnan(val) && first(bounds) ≤ val ≤ last(bounds)
end

function broken_observable_defaults(cf, state=get_defaults_or_inits_dict(cf), defaults=get_defaults_dict(cf))
    obs_defaults = filter(p -> p.first ∈ obssym(cf), pairs(defaults))
    vals = get_initial_state(cf, state, keys(obs_defaults); missing_val=NaN)

    broken = []
    for (val, sym, def) in zip(vals, keys(obs_defaults), values(obs_defaults))
        if !(isapprox(val, def, atol=1e-10))
            push!(broken, (sym, def, val))
        end
    end
    broken
end


initialize_docstring = raw"""
    initialize_componentwise[!](
        nw::Network;
        default_overrides=nothing,
        guess_overrides=nothing,
        bound_overrides=nothing,
        additional_initformula=nothing,
        additional_initconstraint=nothing,
        verbose=false,
        subverbose=false,
        tol=1e-10,
        nwtol=1e-10,
        t=NaN
    ) :: NWState

Initialize a network by solving initialization problems for each component individually,
then verifying the combined solution works for the full network.

There are two version of that function: a mutating one (!-at the end of name) and a non-mutating version.
The mutationg version uses `initialize_component!` internally, the non-mutating one `initialize_component`.
When the mutating version is used, `NWState(nw)` after initialization will return the same initialized
state again, as it is stored in the metadata.

## Parameters
- `nw`: The network to initialize
- `default_overrides`: Dictionary mapping symbolic indices to values that should be used as defaults.
  Use `nothing` as a value for any key to remove that default.
- `guess_overrides`: Dictionary mapping symbolic indices to values to use as initial guesses.
  Use `nothing` as a value for any key to remove that guess.
- `bound_overrides`: Dictionary mapping symbolic indices to bounds for constrained variables.
  Use `nothing` as a value for any key to remove those bounds.
- `additional_initformula`: Dictionary mapping component indices (VIndex/EIndex) to additional initialization formulas.
- `additional_initconstraint`: Dictionary mapping component indices (VIndex/EIndex) to additional initialization constraints.
- `verbose`: Whether to print information about each component initialization
- `subverbose`: Whether to print detailed information within component initialization. Can be Vector [VIndex(1), EIndex(3), ...] for selective output
- `tol`: Tolerance for individual component residuals
- `nwtol`: Tolerance for the full network residual
- `t`: Time at which to evaluate the system

## Returns
- `NWState`: A fully initialized network state that can be used for simulation

## Example of two-step initialization
```julia
# First solve a static model
static_model = create_static_network(...)
static_state = find_fixpoint(static_model)

# Extract interface values and use them to initialize dynamic model
interface_vals = interface_values(static_state)
dynamic_model = create_dynamic_network(...)
dyn_state = initialize_componentwise(dynamic_model, default_overrides=interface_vals)

# Simulate the dynamic model from this initialized state
prob = ODEProblem(dynamic_model, uflat(dyn_state), tspan, pflat(dyn_state))
sol = solve(prob)
```

See also: [`initialize_component`](@ref), [`interface_values`](@ref), [`find_fixpoint`](@ref)
"""
@doc initialize_docstring
initialize_componentwise!(nw; kwargs...) = _initialize_componentwise(Val(true), nw; kwargs...)
@doc initialize_docstring
initialize_componentwise(nw; kwargs...) = _initialize_componentwise(Val(false), nw; kwargs...)
function _initialize_componentwise(
    mutating,
    nw::Network;
    default_overrides=nothing,
    guess_overrides=nothing,
    bound_overrides=nothing,
    additional_initformula=nothing,
    additional_initconstraint=nothing,
    verbose=false,
    subverbose=false,
    tol=1e-10,
    nwtol=1e-10,
    t=NaN
)
    initfun = if mutating == Val(true)
        initialize_component!
    else
        initialize_component
    end
    fullstate = if mutating == Val(true)
        nothing
    else
        Dict{SymbolicIndex, Float64}()
    end
    for vi in 1:nv(nw)
        _default_overrides = _filter_overrides(nw, VIndex(vi), default_overrides)
        _guess_overrides = _filter_overrides(nw, VIndex(vi), guess_overrides)
        _bound_overrides = _filter_overrides(nw, VIndex(vi), bound_overrides)
        _subverbose = _determine_subverbose(subverbose, VIndex(vi))
        verbose && println("Initializing vertex $(vi)...")
        substate = initfun(
            nw[VIndex(vi)],
            default_overrides=_default_overrides,
            guess_overrides=_guess_overrides,
            bound_overrides=_bound_overrides,
            additional_initformula=isnothing(additional_initformula) ? nothing : get(additional_initformula, VIndex(vi), nothing),
            additional_initconstraint=isnothing(additional_initconstraint) ? nothing : get(additional_initconstraint, VIndex(vi), nothing),
            verbose=_subverbose,
            t=t,
            tol=tol
        )
        verbose && _subverbose && println()
        _merge_wrapped!(fullstate, substate, VIndex(vi))
    end
    for ei in 1:ne(nw)
        _default_overrides = _filter_overrides(nw, EIndex(ei), default_overrides)
        _guess_overrides = _filter_overrides(nw, EIndex(ei), guess_overrides)
        _bound_overrides = _filter_overrides(nw, EIndex(ei), bound_overrides)
        _subverbose = _determine_subverbose(subverbose, EIndex(ei))
        verbose && println("Initializing edge $(ei)...")
        substate = initfun(
            nw[EIndex(ei)],
            default_overrides=_default_overrides,
            guess_overrides=_guess_overrides,
            bound_overrides=_bound_overrides,
            additional_initformula=isnothing(additional_initformula) ? nothing : get(additional_initformula, EIndex(ei), nothing),
            additional_initconstraint=isnothing(additional_initconstraint) ? nothing : get(additional_initconstraint, EIndex(ei), nothing),
            verbose=_subverbose,
            t=t,
            tol=tol
        )
        verbose && _subverbose && ei != ne(nw) && println()
        _merge_wrapped!(fullstate, substate, EIndex(ei))
    end

    s0 = if mutating == Val(true)
        NWState(nw)
    else
        usyms = SII.variable_symbols(nw)
        psyms = map(_XPIndex_to_XIndex, SII.parameter_symbols(nw))
        _uflat = [fullstate[s] for s in usyms]
        _pflat = [fullstate[s] for s in psyms]
        NWState(nw, _uflat, _pflat)
    end

    # Calculate the residual for the initialized state
    du = ones(length(uflat(s0)))
    nw(du, uflat(s0), pflat(s0), t)
    resid = LinearAlgebra.norm(du)

    if !(resid < nwtol)
        error("Initialized network has a residual larger than $nwtol: $(resid)! \
               Fix initialization or increase tolerance to suppress error.")
    end
    verbose && println("Initialized network with residual $(resid)!")
    s0
end
_XPIndex_to_XIndex(idx::VPIndex) = VIndex(idx.compidx, idx.subidx)
_XPIndex_to_XIndex(idx::EPIndex) = EIndex(idx.compidx, idx.subidx)
_filter_overrides(_, _, ::Nothing) = nothing
function _filter_overrides(nw, filteridx::SymbolicIndex{Int,Nothing}, dict::AbstractDict)
    filtered = Dict{Symbol, valtype(dict)}()
    for (key, val) in dict
        if filteridx == _baseT(key)(resolvecompidx(nw, key))
            if !(key.subidx isa Symbol)
                error("Overwrites musst be provided as SymbolicIndex{...,Symbol}! Got $key instead.")
            end
            filtered[key.subidx] = val
        end
    end
    filtered
end
function _merge_wrapped!(fullstate, substate, wrapper)
    idxconstructor = _baseT(wrapper)
    for (key, val) in substate
        fullstate[idxconstructor(wrapper.compidx, key)] = val
    end
    fullstate
end
_merge_wrapped!(::Nothing, _, _) = nothing
_determine_subverbose(subverbose::Bool, _) = subverbose
_determine_subverbose(subverbose::AbstractVector, idx) = idx ∈ subverbose
_determine_subverbose(subverbose, idx) = idx == subverbose

"""
    interface_values(s::NWState) :: OrderedDict{SymbolicIndex, Float64}

Extract all interface values (inputs and outputs) from a network state and return them as
a dictionary mapping symbolic indices to their values.

This function is particularly useful in two-step initialization workflows where you want to:
1. Solve a simplified static model first (using [`find_fixpoint`](@ref))
2. Use the resulting interface values to initialize a more complex dynamic model componentwise.

In that scenario, use `interface_values` to for the `default_overrides` argument of
[`initialize_componentwise`](@ref).

See also: [`initialize_componentwise`](@ref), [`find_fixpoint`](@ref) and [`initialize_component`](@ref).
"""
function interface_values(s::NWState)
    nw = extract_nw(s)
    interface_syms = SymbolicIndex[]
    for vi in 1:nv(nw)
        cm = nw[VIndex(vi)]
        syms  = VIndex.(vi, Iterators.flatten((insym_flat(cm), outsym_flat(cm))))
        append!(interface_syms, syms)
    end
    for ei in 1:ne(nw)
        cm = nw[EIndex(ei)]
        syms = EIndex.(ei, Iterators.flatten((insym_flat(cm), outsym_flat(cm))))
        append!(interface_syms, syms)
    end
    @assert allunique(interface_syms) "Interface symbols should be unique!"
    OrderedDict(interface_syms .=> s[interface_syms])
end

"""
    set_interface_defaults!(nw::Network, s::NWState; verbose=false)

Sets the **interface** (i.e., node and edge inputs/outputs) defaults of a given
network to the ones defined by the given state. Notably, while the graph
topology and interface dimensions of the target network `nw` and the source
network of `s` must be identical, the systems may differ in the dynamical
components.

This is mainly intended for initialization purposes: solve the interface values
with a simpler -- possibly static -- network and "transfer" the steady state
interface values to the full network.
"""
function set_interface_defaults!(nw::Network, s::NWState; verbose=false)
    @argcheck s.nw.im.g == nw.im.g "Graphs musst have the same structure!"
    verbose && println("Setting the interface defaults:")
    values = Dict{SymbolicIndex, Float64}()
    for (idxs, _Index, comp) in ((1:nv(nw), VIndex, "Vertex"), (1:ne(nw), EIndex, "Edge"))
        for i in idxs
            target = nw[_Index(i)]
            source = s.nw[_Index(i)]

            tsyms = _Index(i, outsym_flat(target)) |> collect
            svals = s[_Index(i, outsym_flat(source))]
            @assert length(tsyms) == length(svals) "$comp $i: Source network has $(length(svals)) outputs, but target network has $(length(tsyms)) outputs!"
            for (tsym, sval) in zip(tsyms, svals)
                values[tsym] = sval
                verbose && println("  $comp $(lpad(i, length(repr(idxs[end])))) \
                                    output: $(tsym.subidx) => $sval")
            end

            tsyms = _Index(i, insym_flat(target)) |> collect
            svals = s[_Index(i, insym_flat(source))]
            @assert length(tsyms) == length(svals) "$comp $i: Source network has $(length(svals)) inputs, but target network has $(length(tsyms)) inputs!"
            for (tsym, sval) in zip(tsyms, svals)
                values[tsym] = sval
                verbose && println("  $comp $(lpad(i, length(repr(idxs[end])))) \
                                    input:  $(tsym.subidx) => $sval")
            end
        end
    end
    for (sym, val) in values
        set_default!(nw, sym, val)
    end
    nw
end
