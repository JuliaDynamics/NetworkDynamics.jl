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

function initialization_problem(cf::T;
    t=NaN ,
    defaults=get_defaults_dict(cf),
    guesses=get_guesses_dict(cf),
    bounds=get_bounds_dict(cf),
    apply_bound_transformation,
    verbose=true) where {T<:ComponentModel}
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

    Neqs = dim(cf)  + mapreduce(length, +, outfree_ms)

    freesym = vcat(mapreduce((syms, map) -> syms[map], vcat, outsym_normalized(cf), outfree_ms),
                   sym(cf)[ufree_m],
                   mapreduce((syms, map) -> syms[map], vcat, insym_normalized(cf), infree_ms),
                   psym(cf)[pfree_m])

    @assert length(freesym) == Nfree
    if Neqs < Nfree
        throw(ArgumentError("Initialization problem underconstraint. $(Neqs) Equations for $(Nfree) free variables: $freesym"))
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

    N = ForwardDiff.pickchunksize(Nfree)
    fg = compfg(cf)
    unlcache = map(d->DiffCache(zeros(d), N), length(freesym))
    outcaches = map(d->DiffCache(zeros(d), N), outdim_normalized(cf))
    ucache = DiffCache(zeros(dim(cf)), N)
    incaches = map(d->DiffCache(zeros(d), N), indim_normalized(cf))
    pcache = DiffCache(zeros(pdim(cf)))

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
        @views dunl_out = dunl[dim(cf)+1:end]

        # this fills the second half of the du buffer with the fixed and current outputs
        dunl_out .= RecursiveArrayTools.ArrayPartition(outbufs...)
        # execute fg to fill dunl and outputs
        fg(outbufs, dunl_fg, ubuf, inbufs, pbuf, t)

        # calculate the residual for the second half ov the dunl buf, the outputs
        dunl_out .= dunl_out .- RecursiveArrayTools.ArrayPartition(outbufs...)
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
                         additional_defaults=nothing,
                         additional_guesses=nothing,
                         additional_bounds=nothing,
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
- `additional_defaults/guesses/bounds`: Dictionary to merge with `defaults`/`guesses`/`bounds`
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
                             additional_defaults=nothing,
                             additional_guesses=nothing,
                             additional_bounds=nothing,
                             verbose=true,
                             apply_bound_transformation=true,
                             t=NaN,
                             tol=1e-10,
                             kwargs...)

    defaults = isnothing(additional_defaults) ? defaults : merge(defaults, additional_defaults)
    guesses  = isnothing(additional_guesses)  ? guesses  : merge(guesses, additional_guesses)
    bounds   = isnothing(additional_bounds)   ? bounds   : merge(bounds, additional_bounds)

    prob, boundT! = initialization_problem(cf; defaults, guesses, bounds, verbose, apply_bound_transformation)

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
        verbose && @info "No free variables! Residual $(LinearAlgebra.norm(residual))"
    end
    if res > tol
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
                          verbose=true,
                          t=NaN,
                          kwargs...)

Mutating version of [`initialize_component`](@ref). See this docstring for all details.
In contrast to the non mutating version, this function reads in defaults and guesses from the
symbolic metadata and writes the initialized values back in to the metadata.

## Parameters
- `cf`: ComponentModel to initialize
- `defaults`: Optional dictionary to replace metadata defaults
- `guesses`: Optional dictionary to replace metadata guesses
- `bounds`: Optional dictionary to replace metadata bounds
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
                              verbose=true,
                              t=NaN,
                              kwargs...)

    # Synchronize defaults, guesses, and bounds
    _sync_metadata!(cf, defaults, get_defaults_dict, set_default!, delete_default!, "default", verbose)
    _sync_metadata!(cf, guesses,  get_guesses_dict,  set_guess!,   delete_guess!,   "guess",   verbose)
    _sync_metadata!(cf, bounds,   get_bounds_dict,   set_bounds!,  delete_bounds!,  "bounds",  verbose)

    # Now proceed with initialization using the updated metadata
    init_state = initialize_component(
        cf;
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
        has_default(cf, sym) && continue
        set_init!(cf, sym, val)
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
        sym(cf),                               # states
        filter(s->is_unused(cf, s), psym(cf)), # used parameters
        insym_flat(cf),                        # inputs
        outsym_flat(cf)                        # outputs
    ))

    if keys(state) ⊆ needed_symbols
        missing = setdiff(needed_symbols, keys(state))
        throw(ArgumentError("State dictionary is missing required symbols: $missing. \
                             These symbols need values to calculate the residual."))
    end

    outs = Tuple(Float64[get(state, s, NaN) for s in sv] for sv in outsym_normalized(cf))
    u = Float64[get(state, s, NaN) for s in sym(cf)]
    ins = Tuple(Float64[get(state, s, NaN) for s in sv] for sv in insym_normalized(cf))
    p = Float64[is_unused(cf, s) ? NaN : get(state, s, NaN) for s in psym(cf)]

    res = zeros(dim(cf) + sum(outdim_normalized(cf)))

    res_fg = @views res[1:dim(cf)]
    res_out = @views res[dim(cf)+1:end]
    res_out .= RecursiveArrayTools.ArrayPartition(outs...)
    compfg(cf)(outs, res_fg, u, ins, p, t)
    res_out .= res_out .- RecursiveArrayTools.ArrayPartition(outs...)

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
        if !(val ≈ def)
            push!(broken, (sym, def, val))
        end
    end
    broken
end

"""
    set_interface_defaults!(nw::Network, s::NWState; verbose=false)

Sets the **interface** (i.e. node and edge inputs/outputs) defaults of a given
network to the ones defined by the given state. Notably, while the graph
topology and interface dimensions of the target network `nw` and the source
network of `s` musst be identicaly, the systems may differ in the dynamical
components.

This is mainly inteded for initialization purposes: solve the interface values
with a simpler -- possible static -- network and "transfer" the steady state
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
