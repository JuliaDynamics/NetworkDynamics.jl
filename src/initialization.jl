struct NetworkInitError <: Exception
    msg::String
end
struct ComponentInitError <: Exception
    msg::String
end
Base.showerror(io::IO, e::NetworkInitError) = print(io, "NetworkInitError: ", e.msg)
Base.showerror(io::IO, e::ComponentInitError) = print(io, "ComponentInitError: ", e.msg)

"""
    find_fixpoint(nw::Network, [x0::NWState=NWState(nw)], [p::NWParameter=x0.p]; kwargs...)
    find_fixpoint(nw::Network, x0::AbstractVector, p::AbstractVector; kwargs...)

Find a steady-state (fixed-point) solution of the network dynamics by solving the nonlinear
equation `f(u, p, t) = 0`, where `f` represents the network's right-hand side function.

This is a convenience wrapper around `SteadyStateProblem` from the SciML ecosystem that
constructs and solves the steady state problem, returning the solution as an `NWState`.

## Arguments
- `nw::Network`: The network dynamics to find a fixed point for
- `x0`: Initial guess for the state variables. Can be:
  - `NWState`: Complete network state (default: `NWState(nw; ufill=0)`)
  - `AbstractVector`: Flat state vector
- `p`: Network parameters. Can be:
  - `NWParameter`: Complete parameter object (default: extracted from `x0` or `NWParameter(nw)`)
  - `AbstractVector`: Flat parameter vector

## Keyword Arguments
- `alg=SSRootfind()`: Steady state solver algorithm from NonlinearSolve.jl
- `t=NaN`: Time at which to evaluate the system (for time-dependent networks)
- Additional `kwargs` are passed to the SciML `solve` function

## Returns
- `NWState`: Network state at the found fixed point

See also: [`NWState`](@ref), [`NWParameter`](@ref), [`initialize_componentwise`](@ref)
"""
function find_fixpoint(nw::Network, x0::AbstractVector; kwargs...)
    find_fixpoint(nw, NWState(nw, x0), NWParameter(nw); kwargs...)
end
function find_fixpoint(nw::Network, p::NWParameter; kwargs...)
    find_fixpoint(nw, NWState(p; ufill=0), p; kwargs...)
end
function find_fixpoint(nw::Network,
                       x0::NWState=NWState(nw; ufill=0),
                       p::NWParameter=x0.p;
                       t=isnothing(x0.t) ? NaN : x0.t,
                       kwargs...)
    find_fixpoint(nw, uflat(x0), pflat(p); t, kwargs...)
end

#=
We want to fix time for steady state to allow for initializing at varios time points (SteadyStateDiffEq uses t->inf)
However, we build a custom struct for that to still enable symbolic indexing and sparsity!
=#
struct NetworkFixedT{NW}
    nw::NW
    t::Float64
    NetworkFixedT(nw, t) = new{typeof(nw)}(nw, Float64(t))
end
(nwft::NetworkFixedT)(du, u, p, _) = nwft.nw(du, u, p, nwft.t)
SciMLBase.__has_sys(nw::NetworkFixedT) = true
SciMLBase.__has_jac_prototype(nw::NetworkFixedT) = !isnothing(nw.nw.jac_prototype)
function Base.getproperty(nw::NetworkFixedT, s::Symbol)
    if s===:sys
        nw.nw
    elseif s===:jac_prototype
        getfield(nw.nw, :jac_prototype)[]
    else
        getfield(nw, s)
    end
end

function find_fixpoint(nw::Network, x0::AbstractVector, p::AbstractVector;
                       alg=SSRootfind(), t=NaN, kwargs...)
    nw_fixed_t = NetworkFixedT(nw, t)

    du = zeros(length(x0))
    fn_error = nothing
    try
        nw_fixed_t(du, x0, p, nothing)
    catch e
        fn_error = e
    end
    output_nans = any(isnan, du)

    if !isnothing(fn_error) || output_nans
        # something bad happened
        io = IOBuffer()
        if isnothing(fn_error)
            nans = findall(isnan, du)
            syms = SII.variable_symbols(nw)[nans]
            println(io, "Evaluation of network rhs at t=$t resulted in NaNs.)")
            println(io, "  NaN-Results:  ", join(syms, ", "))
            print(io,"Possible causes:")
        else
            print(io, "Evaluation of network rhs at t=$t led to an error! Possible causes:")
        end

        if any(isnan, x0) || any(isnan, p)
            print(io, "\n - Some input states/parameters contained NaN:")
            if any(isnan, x0)
                nan_indices = findall(isnan, x0)
                variable_symbols = NetworkDynamics.SII.variable_symbols(nw)
                missing_vars = [string(variable_symbols[i]) for i in nan_indices]
                print(io, "\n   Variables:  ", join(missing_vars, ", "))
            end

            if any(isnan, p)
                nan_indices = findall(isnan, p)
                parameter_symbols = NetworkDynamics.SII.parameter_symbols(nw)
                missing_params = [string(parameter_symbols[i]) for i in nan_indices]
                print(io, "\n   Parameters: ", join(missing_params, ", "))
            end
        end

        if output_nans && isnan(t)
            print(io, "\n - Network was evaluated at t=NaN, maybe your system has explicit time dependence? Try specifying kw argument `t` to decide on time")
        end
        if !isnothing(fn_error)
            print(io, "\nCausing ERROR:\n")
            Base.showerror(io, fn_error)
        end

        throw(NetworkInitError(String(take!(io))))
    end

    prob = SteadyStateProblem(nw_fixed_t, x0, p)
    sol = SciMLBase.solve(prob, alg; kwargs...)
    if !SciMLBase.successful_retcode(sol.retcode)
        throw(NetworkInitError("""
        Could not find fixpoint, solver returned $(sol.retcode) (alg=$(alg))! For debugging, \
        it is advised to manually construct the steady state problem and try \
        different solvers/arguments:

        prob = SteadyStateProblem(nw, uflat(nwstate), pflat(nwpara))
        sol = solve(prob, alg; kwargs...)
        x0 = NWState(nw, sol.u)

        For detail see https://docs.sciml.ai/NonlinearSolve/stable/native/steadystatediffeq/
        """))
    end
    NWState(nw, sol.u, NWParameter(nw, p))
end

function initialization_problem(cf::T,
    defaults,
    guesses,
    bounds,
    initconstraint=nothing;
    t=NaN,
    apply_bound_transformation=true,
    verbose=true,
    io=stdout
) where {T<:ComponentModel}
    hasinsym(cf) || throw(ArgumentError("Component model must have `insym` set to support initialiation!"))

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
            printstyled(io, " - Apply positivity/negativity conserving variable transformation on $(freesym[idxs]) to satisfy bounds.\n")
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
    isempty(missing_guesses) || throw(ComponentInitError("Missing guesses for free variables $(missing_guesses)"))

    # apply bound conserving transformation to initial state
    try
        inv_boundT!(uguess)
    catch e
        throw(ComponentInitError("Failed to apply bound-conserving transformation to initial guess. \
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
        if verbose
            printstyled(io, " - Initialization problem is fully constrained. Created NonlinearLeastSquaresProblem for:\n")
            for (sym, guess) in zip(freesym, uguess)
                print(io, "   - ", sym, " (guess=$(guess))\n")
            end
        end
    elseif Neqs > Nfree
        if verbose
            printstyled(io, " - Initialization problem is overconstrained ($Nfree vars for $Neqs equations). Create NonlinearLeastSquaresProblem for:\n")
            for (sym, guess) in zip(freesym, uguess)
                print(io, "   - ", sym, " (guess=$(guess))\n")
            end
        end
    else
        # verbose && printstyled(io, "Initialization problem is underconstrained ($Nfree vars for $Neqs equations). Create NonlinearLeastSquaresProblem for $freesym.\n"; color=:yellow)
        printstyled(io, " - WARNING:", color=:yellow)
        printstyled(io, "Initialization problem is underconstrained ($Nfree vars for $Neqs equations). Create NonlinearLeastSquaresProblem for:")
        for (sym, guess) in zip(freesym, uguess)
            print(io, "   - ", sym, " (guess=$(guess))\n")
        end
    end

    # check rhs of system for obvious problems
    resid = zeros(Neqs)
    fn_error = nothing
    try
        nlf(resid, uguess, nothing)
    catch e
        fn_error = e
    end
    resid_nan = any(isnan, resid)
    if !isnothing(fn_error) || resid_nan
        io = IOBuffer()
        print(io, "Error while constructing initialization problem for $(cf.name):")
        if resid_nan
            print(io, " - Residual contains NaNs!\n")
            if isnan(t)
                print(io, " System initialized at t=NaN, maybe your system has explicit \
                    time dependence? Try specifying kw argument `t` to decide on time")
            end
        end
        if !isnothing(fn_error)
            print(io, " - Error while calling RHS of initialization problem!\n")
            print(io, "\n\nOriginal Error:\n")
            Base.showerror(io, fn_error)
        end
        throw(ComponentInitError(String(take!(io))))
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
                         additional_guessformula=nothing,
                         additional_initconstraint=nothing,
                         verbose=true,
                         apply_bound_transformation=true,
                         t=NaN,
                         tol=1e-10,
                         residual=nothing,
                         alg=nothing, # defaults to FastShortcutNLLSPolyalg(linsolve=QRFactorization())
                         solve_kwargs=(;),
                         io=stdout,
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
- `additional_initformula`: Additional initialization formulas to apply beyond those in component metadata.
  InitFormulas compute and set default values, reducing the number of free variables.
- `additional_guessformula`: Additional guess formulas to apply beyond those in component metadata.
  GuessFormulas compute improved initial guesses for free variables, improving solver convergence.
- `additional_initconstraint`: Additional initialization constraints to apply beyond those in component metadata
- `verbose`: Whether to print information during initialization
- `apply_bound_transformation`: Whether to apply bound-conserving transformations
- `t`: Time at which to solve for steady state. Only relevant for components with explicit time dependency.
- `tol`: Tolerance for the residual of the initialized model (defaults to `1e-10`). Init throws error if resid ≥ tol.
- `residual`: Optional `Ref{Float64}` which gets the final residual of the initialized model.
- `alg=nothing`: Nonlinear solver algorithm (defaults to NonlinearSolve.jl default with QR factorization, since init problems tend to be ill-conditioned.)
- `solve_kwargs=(;)`: Additional keyword arguments passed to the SciML `solve` function
- `io=stdout`: IO stream for printing information

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
                             additional_guessformula=nothing,
                             additional_initconstraint=nothing,
                             verbose=true,
                             apply_bound_transformation=true,
                             t=NaN,
                             tol=1e-10,
                             residual=nothing,
                             # internal keywords to "return" the final defaults/guesses after applying formulas
                             _final_defaults=nothing,
                             _final_guesses=nothing,
                             alg=nothing,
                             solve_kwargs=(;),
                             io=stdout,
                             kwargs...)

    if !isempty(kwargs)
        @warn "Passing `kwargs` to `initialize_component(!)` is deprecated. Use `alg` and `solve_kwargs=(; kw=val)` instead."
    end

    if alg == nothing
        alg = FastShortcutNLLSPolyalg(linsolve=QRFactorization())
    end

    defaults = isnothing(default_overrides) ? defaults : merge(defaults, default_overrides)
    guesses  = isnothing(guess_overrides)   ? guesses  : merge(guesses, guess_overrides)
    bounds   = isnothing(bound_overrides)   ? bounds   : merge(bounds, bound_overrides)

    # filter out nothing values
    defaults = filter(p -> !isnothing(p.second), defaults)
    guesses = filter(p -> !isnothing(p.second), guesses)
    bounds = filter(p -> !isnothing(p.second), bounds)

    # Extract metadata and merge with additional constraints/formulas
    metadata_initformulas = has_initformula(cf) ? get_initformulas(cf) : nothing
    combined_initformulas = collect_formulas(metadata_initformulas, additional_initformula)

    # Apply initialization formulas to defaults
    if !isnothing(combined_initformulas)
        apply_init_formulas!(defaults, combined_initformulas; verbose, io)
    end

    # Extract and apply guess formulas
    metadata_guessformulas = has_guessformula(cf) ? get_guessformulas(cf) : nothing
    combined_guessformulas = collect_formulas(metadata_guessformulas, additional_guessformula)
    if !isnothing(combined_guessformulas)
        apply_guess_formulas!(guesses, defaults, combined_guessformulas; verbose, io)
    end

    metadata_constraint = has_initconstraint(cf) ? get_initconstraints(cf) : nothing
    combined_constraint = merge_initconstraints(metadata_constraint, additional_initconstraint)

    # optinally "return" the final defaults and guesses via keyword ref
    if _final_defaults isa Ref
        _final_defaults[] = defaults
    end
    if _final_guesses isa Ref
        _final_guesses[] = guesses
    end

    prob, boundT! = initialization_problem(cf, defaults, guesses, bounds,
                                          combined_constraint;
                                          verbose, apply_bound_transformation, t, io)

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
        sol = SciMLBase.solve(prob, alg; verbose, kwargs..., solve_kwargs...)

        if sol.prob isa NonlinearLeastSquaresProblem && sol.retcode == SciMLBase.ReturnCode.Stalled
            res = LinearAlgebra.norm(sol.resid)
            printstlyled(" - WARN: "; color=:yellow)
            printstyled("Initialization for component stalled with residual $(res)")
        elseif !SciMLBase.successful_retcode(sol.retcode)
            throw(ComponentInitError("Initialization failed. Solver returned $(sol.retcode)"))
        else
            res = LinearAlgebra.norm(sol.resid)
            verbose && printstyled(io, " - Initialization successful with residual $(res)\n")
        end

        # Transform back to original space
        u = boundT!(copy(sol.u))

        # Add solved values to the complete dictionary
        free_symbols = SII.variable_symbols(sol)
        for (sym, val) in zip(free_symbols, u)
            verbose && print(io, "   - ", sym, " => ", val,"\n")
            init_state[sym] = val
        end
    else
        res = init_residual(cf, init_state; t)
        verbose && printstyled(io, " - No free variables! Residual $(res)\n")
    end
    if residual isa Ref
        residual[] = res
    end
    if !(res < tol)
        throw(ComponentInitError("Initialized model has a residual larger than specified tolerance $(res) > $(tol)! \
               Fix initialization or increase tolerance to suppress error."))
    end


    # Check for broken bounds using the complete state
    broken_bnds = broken_bounds(cf, init_state, bounds)
    if !isempty(broken_bnds)
        broken_msgs = ["  $sym = $val (bounds: $lb..$ub)" for (sym, val, (lb, ub)) in broken_bnds]
        fullmsg = "Initialized model has broken bounds. Try to adapt the initial guesses!" *
              "\n" * join(broken_msgs, "\n")
        if verbose
            printstyled(io, " - WARNING: ", color=:yellow)
            printstyled(io, fullmsg * "\n")
        else
            @warn fullmsg
        end
    end

    # Check for broken observable defaults
    if !isempty(observable_defaults)
        broken_obs = broken_observable_defaults(cf, init_state, observable_defaults)
        if !isempty(broken_obs)
            broken_msgs = ["  $sym = $val (default: $def)" for (sym, def, val) in broken_obs]
            fullmsg = "Initialized model has observables that differ from their specified defaults:" *
                  "\n" * join(broken_msgs, "\n")
            if verbose
                printstyled(io, " - WARNING: ", color=:yellow)
                printstyled(io, fullmsg * "\n")
            else
                @warn fullmsg
            end
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
                          additional_guessformula=nothing,
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
- `additional_initformula`: Additional initialization formulas to apply beyond those in component metadata.
  InitFormulas compute and set default values, reducing the number of free variables.
- `additional_guessformula`: Additional guess formulas to apply beyond those in component metadata.
  GuessFormulas compute improved initial guesses for free variables, improving solver convergence.
- `additional_initconstraint`: Additional initialization constraints to apply beyond those in component metadata
- `verbose`: Whether to print information during initialization
- `t`: Time at which to solve for steady state. Only relevant for components with explicit time dependency.
- All other `kwargs` are passed to [`initialize_component`](@ref), see its docstring for details!

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
                              additional_guessformula=nothing,
                              additional_initconstraint=nothing,
                              verbose=true,
                              t=NaN,
                              io=stdout,
                              kwargs...)

    # Synchronize defaults, guesses, and bounds
    _sync_metadata!(cf, defaults, get_defaults_dict, set_default!, delete_default!, "default", verbose, io)
    _sync_metadata!(cf, guesses,  get_guesses_dict,  set_guess!,   delete_guess!,   "guess",   verbose, io)
    _sync_metadata!(cf, bounds,   get_bounds_dict,   set_bounds!,  delete_bounds!,  "bounds",  verbose, io)

    for (name, set_fn!, rm_fn!, dict) in (
        ("default", set_default!, delete_default!, default_overrides),
        ("guess", set_guess!, delete_guess!, guess_overrides),
        ("bound", set_bounds!, delete_bounds!, bound_overrides)
    )
        isnothing(dict) && continue
        for (sym, val) in dict
            if !isnothing(val)
                set_fn!(cf, sym, val)
                verbose && println(io, " - Set additional $name for $sym: $val")
            else
                rm_fn!(cf, sym)
                verbose && println(io, " - Remove $name for $sym")
            end
        end
    end

    # make sure to get back the final defaults/guesses after applying formulas
    _final_defaults = Ref{Dict{Symbol,Float64}}()
    _final_guesses = Ref{Dict{Symbol,Float64}}()

    # Now proceed with initialization using the updated metadata
    init_state = initialize_component(
        cf;
        additional_initformula=additional_initformula,
        additional_guessformula=additional_guessformula,
        additional_initconstraint=additional_initconstraint,
        verbose=verbose,
        t=t,
        _final_defaults,
        _final_guesses,
        io,
        kwargs...  # Only pass the remaining kwargs
    )

    # delete all inits
    for s in keys(get_inits_dict(cf))
        delete_init!(cf, s)
    end

    # write back defaults/guesses if they've changed due to formulas
    for (sym, val) in _final_defaults[]
        if !has_default(cf, sym) || get_default(cf, sym) != val
           set_default!(cf, sym, val)
        end
    end
    for (sym, val) in _final_guesses[]
        if !has_guess(cf, sym) || get_guess(cf, sym) != val
            set_guess!(cf, sym, val)
        end
    end
    # Write back the initialized values to the component metadata
    for (sym, val) in init_state
        if !has_default(cf, sym)
            # No default exists, store as init value
            set_init!(cf, sym, val)
        end
    end

    cf
end
function _sync_metadata!(cf, provided, get_orig_fn, set_fn!, remove_fn!, type_name, verbose, io)
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
        verbose && println(io, "Removed $type_name for $sym (was $(original[sym]))")
    end

    # Add new keys
    for sym in new_keys
        set_fn!(cf, sym, provided[sym])
        verbose && println(io, "Added $type_name for $sym: $(provided[sym])")
    end

    # Update changed keys
    for sym in changed_keys
        set_fn!(cf, sym, provided[sym])
        verbose && println(io, "Updated $type_name for $sym: $(original[sym]) → $(provided[sym])")
    end
end


function isinitialized(cf::ComponentModel)
    all(has_default_or_init(cf, s) || is_unused(cf, s) for s in vcat(sym(cf), psym(cf)))
end

"""
    init_residual(cf::ComponentModel, [state=get_defaults_or_inits_dict(cf)]; t=NaN, verbose=false)

Calculates the residual |du| for the given component model using the values provided.
If no state dictionary is provided, it uses the values from default/init [Metadata](@ref).

If `verbose=true` prints the residual of every single state.

See also [`initialize_component`](@ref).
"""
function init_residual(cf::ComponentModel, state=get_defaults_or_inits_dict(cf); t=NaN, verbose=false)
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
    fill!(res, NaN)

    res_fg = @views res[1:dim(cf)]
    res_out = @views res[dim(cf)+1:dim(cf)+Nout]
    res_add = @views res[dim(cf)+Nout+1:end]

    # buffer to store the calculated output at point
    out_calc = map(outs) do out
        buf = similar(out)
        fill!(buf, NaN) # fill with NaN
    end
    compfg(cf)(out_calc, res_fg, u, ins, p, t)

    # difference between calculated outputs and provided outputs
    res_out .= RecursiveArrayTools.ArrayPartition(out_calc...) .- RecursiveArrayTools.ArrayPartition(outs...) # compare to calculated ouputs

    additional_cf(res_add, outs, u, ins, p, t) # apply additional constraints

    if verbose
        _sym = sym(cf)
        _osym = outsym(cf)
        _addsym = ["init constraint $i" for i in 1:additional_Neqs]
        lines = String[]
        for (i, s) in enumerate(vcat(_sym, _osym, _addsym))
            push!(lines, "  $s &=> $(res[i])")
        end
        aligned = align_strings(lines)
        for l in aligned
            println(l)
        end
    end

    res =  LinearAlgebra.norm(res)
    if isnan(res)
        err_str = if isnan(t)
            "Residual of component is NaN at t=NaN! Maybe your system has \
                explicit time dependence? Try specifying kw argument `t` to decide on time."
        else
            "Residual of component is NaN, which should not happen for an initialized system!"
        end
        throw(ComponentInitError(err_str))
    end
    return res
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
        additional_guessformula=nothing,
        additional_initconstraint=nothing,
        verbose=false,
        subverbose=false,
        tol=1e-10,
        nwtol=1e-10,
        t=NaN,
        subalg=nothing,
        subsolve_kwargs=nothing,
        parallel=false,
    ) :: NWState

Initialize a network by solving initialization problems for each component individually,
then verifying the combined solution works for the full network.

There are two versions of that function: a mutating one (!-at the end of name) and a non-mutating version.
The mutating version uses `initialize_component!` internally, the non-mutating one `initialize_component`.
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
  InitFormulas compute and set default values, reducing the number of free variables.
- `additional_guessformula`: Dictionary mapping component indices (VIndex/EIndex) to additional guess formulas.
  GuessFormulas compute improved initial guesses for free variables, improving solver convergence.
- `additional_initconstraint`: Dictionary mapping component indices (VIndex/EIndex) to additional initialization constraints.
- `verbose`: Whether to print information about each component initialization
- `subverbose`: Whether to print detailed information within component initialization. Can be Vector [VIndex(1), EIndex(3), ...] for selective output
- `tol`: Tolerance for individual component residuals
- `nwtol`: Tolerance for the full network residual
- `t`: Time at which to evaluate the system
- `subalg`: Nonlinear solver algorithm to use for component initialization (defaults to NonlinearSolve.jl default). Can be passed as single value or dict mapping VIndex/EIndex to alg (non-existent keys use default).
- `subsolve_kwargs`: Additional keyword arguments passed to the SciML `solve` function for component initialization.
  Can be passed as single value or dict mapping VIndex/EIndex to kwargs (non-existent keys use empty kwargs `(;)`).
- `parallel=false`: (Experimental) Whether to initialize components in parallel using multithreading.

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
initialize_componentwise!(nw; kwargs...) = _initialize_componentwise(Val{true}(), nw; kwargs...)
@doc initialize_docstring
initialize_componentwise(nw; kwargs...) = _initialize_componentwise(Val{false}(), nw; kwargs...)
function _initialize_componentwise(
    mutating,
    nw::Network;
    default_overrides=nothing,
    guess_overrides=nothing,
    bound_overrides=nothing,
    additional_initformula=nothing,
    additional_guessformula=nothing,
    additional_initconstraint=nothing,
    verbose=false,
    subverbose=false,
    tol=1e-10,
    nwtol=1e-10,
    t=NaN,
    subalg=nothing,
    subsolve_kwargs=nothing,
    parallel=false,#(ne(nw) + nv(nw)) > 200,
)
    # Check for aliased components in mutating mode
    if mutating == Val{true}()
        has_aliased_v = !isempty(nw.im.aliased_vertexms)
        has_aliased_e = !isempty(nw.im.aliased_edgems)
        if has_aliased_v || has_aliased_e
            component_type = if has_aliased_v && has_aliased_e
                "vertices and edges"
            elseif has_aliased_v
                "vertices"
            else
                "edges"
            end
            throw(ArgumentError("""
                Cannot use mutating initialization (initialize_componentwise!) with aliased components!

                Your network contains aliased $component_type (the same component object
                is referenced at multiple indices). When using the mutating version,
                initialization would modify the shared component's metadata, causing all
                aliased instances to have identical initialization values.

                Solutions:
                1. Use the non-mutating version: initialize_componentwise(nw, ...)
                2. Reconstruct the network with dealias=true: Network(g, vertexm, edgem; dealias=true)
                3. Manually copy component models before creating the network
                """))
        end
    end

    initfun = if mutating == Val{true}()
        initialize_component!
    else
        initialize_component
    end
    fullstate = if mutating == Val{true}()
        nothing
    else
        Dict{SymbolicIndex, Float64}()
    end
    statelock = ReentrantLock()

    # dict might be provided as VIndex(:foo) so we need to resolve comp names
    subalg = _resolve_subargument_dict(nw, subalg)
    subsolve_kwargs = _resolve_subargument_dict(nw, subsolve_kwargs)

    # to improve type stability, we split the overrides into v/e parts only once
    default_v_overrides, default_e_overrides = _split_overrides(default_overrides)
    guess_v_overrides, guess_e_overrides = _split_overrides(guess_overrides)
    bound_v_overrides, bound_e_overrides = _split_overrides(bound_overrides)

    tasks = SpinTask[]
    for vi in 1:nv(nw)
        task = SpinTask("Initialize Vertex $vi") do io
            _comp = nw[VIndex(vi)]
            _default_overrides = _filter_overrides(nw, VIndex(vi), default_v_overrides)
            _guess_overrides = _filter_overrides(nw, VIndex(vi), guess_v_overrides)
            _bound_overrides = _filter_overrides(nw, VIndex(vi), bound_v_overrides)
            _add_initformula=isnothing(additional_initformula) ? nothing : get(additional_initformula, VIndex(vi), nothing)
            _add_guessformula=isnothing(additional_guessformula) ? nothing : get(additional_guessformula, VIndex(vi), nothing)
            _add_initconstraint=isnothing(additional_initconstraint) ? nothing : get(additional_initconstraint, VIndex(vi), nothing)
            _subverbose = _determine_subverbose(subverbose, VIndex(vi))
            _alg = _determine_subargument(subalg, VIndex(vi), nothing)
            _solve_kwargs = _determine_subargument(subsolve_kwargs, VIndex(vi), (;))
            rescapture = Ref{Float64}(NaN) # residual for the component

            _subverbose && println()
            substate = initfun(
                _comp,
                default_overrides=_default_overrides,
                guess_overrides=_guess_overrides,
                bound_overrides=_bound_overrides,
                additional_initformula=_add_initformula,
                additional_guessformula=_add_guessformula,
                additional_initconstraint=_add_initconstraint,
                verbose=_subverbose,
                t=t,
                tol=tol,
                residual=rescapture,
                alg=_alg,
                solve_kwargs=_solve_kwargs,
                io=io,
            )
            _merge_wrapped!(fullstate, substate, VIndex(vi), statelock)
            # return residual, if subverbose printed anyway
            _subverbose ? nothing : rescapture[]
        end
        push!(tasks, task)
    end
    for ei in 1:ne(nw)
        task = SpinTask("Initialize Edge $ei") do io
            _comp = nw[EIndex(ei)]
            _default_overrides = _filter_overrides(nw, EIndex(ei), default_e_overrides)
            _guess_overrides = _filter_overrides(nw, EIndex(ei), guess_e_overrides)
            _bound_overrides = _filter_overrides(nw, EIndex(ei), bound_e_overrides)
            _add_initformula=isnothing(additional_initformula) ? nothing : get(additional_initformula, EIndex(ei), nothing)
            _add_guessformula=isnothing(additional_guessformula) ? nothing : get(additional_guessformula, EIndex(ei), nothing)
            _add_initconstraint=isnothing(additional_initconstraint) ? nothing : get(additional_initconstraint, EIndex(ei), nothing)
            _subverbose = _determine_subverbose(subverbose, EIndex(ei))
            _alg = _determine_subargument(subalg, EIndex(ei), nothing)
            _solve_kwargs = _determine_subargument(subsolve_kwargs, EIndex(ei), (;))
            rescapture = Ref{Float64}(NaN) # residual for the component

            _subverbose && println()
            substate = initfun(
                _comp,
                default_overrides=_default_overrides,
                guess_overrides=_guess_overrides,
                bound_overrides=_bound_overrides,
                additional_initformula=_add_initformula,
                additional_guessformula=_add_guessformula,
                additional_initconstraint=_add_initconstraint,
                verbose=_subverbose,
                t=t,
                tol=tol,
                residual=rescapture,
                alg=_alg,
                solve_kwargs=_solve_kwargs,
                io=io,
            )
            _merge_wrapped!(fullstate, substate, EIndex(ei), statelock)
            # return residual, if subverbose printed anyway
            _subverbose ? nothing : rescapture[]
        end
        push!(tasks, task)
    end

    if parallel
        bthreads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        if stdout isa Base.TTY # check if print backend supports deletion
            run_fancy(tasks; verbose)
        else
            run_plain(tasks; verbose)
        end
        BLAS.set_num_threads(bthreads)
    else
        run_sequential(tasks; verbose)
    end

    s0 = if mutating == Val{true}()
        NWState(nw)
    else
        usyms = SII.variable_symbols(nw)
        psyms = SII.parameter_symbols(nw)
        _uflat = [fullstate[s] for s in usyms]
        _pflat = [fullstate[s] for s in psyms]
        NWState(nw, _uflat, _pflat)
    end

    # Calculate the residual for the initialized state
    du = ones(length(uflat(s0)))
    nw(du, uflat(s0), pflat(s0), t)
    resid = LinearAlgebra.norm(du)

    if !(resid < nwtol)
        throw(NetworkInitError("Initialized network has a residual larger than $nwtol: $(resid)! \
               Fix initialization or increase tolerance to suppress error."))
    end
    verbose && println("Initialized network with residual $(resid)!")
    s0
end
_filter_overrides(_, _, ::Nothing) = nothing
function _filter_overrides(nw, filteridx, dict::AbstractDict)
    filtered = Dict{Symbol, valtype(dict)}()
    for (key, val) in dict
        if filteridx.compidx == resolvecompidx(nw, key)
            filtered[key.subidx] = val
        end
    end
    filtered
end
_split_overrides(::Nothing) = nothing, nothing
function _split_overrides(dict)
    voverrides = Dict{VIndex{Int,Symbol}, valtype(dict)}()
    eoverrides = Dict{EIndex{Int,Symbol}, valtype(dict)}()
    for (key, val) in dict
        if key isa VIndex{Int,Symbol}
            voverrides[key] = val
        elseif key isa EIndex{Int,Symbol}
            eoverrides[key] = val
        else
            error("Overrides must be provided as SymbolicIndex{Int,Symbol} (VIndex or EIndex)! Got $key instead.")
        end
    end
    voverrides, eoverrides
end
function _merge_wrapped!(fullstate, substate, wrapper, lock)
    idxconstructor = idxtype(wrapper)
    @lock lock begin
        for (key, val) in substate
            fullstate[idxconstructor(wrapper.compidx, key)] = val
        end
    end
    nothing
end
_merge_wrapped!(::Nothing, _, _, _) = nothing
_determine_subverbose(subverbose::Bool, _) = subverbose
_determine_subverbose(subverbose::AbstractVector, idx) = idx ∈ subverbose
_determine_subverbose(subverbose, idx) = idx == subverbose

_determine_subargument(subargument::Nothing, idx, default) = default
_determine_subargument(subargument::AbstractDict, idx, default) = get(subargument, idx, default)
_determine_subargument(subargument, idx, default) = subargument

_resolve_subargument_dict(nw, nodict) = nodict # no dict
_resolve_subargument_dict(nw, gooddict::AbstractDict{<:SymbolicIndex{Int}}) = gooddict
function _resolve_subargument_dict(nw, dict::AbstractDict)
    if any(t -> !(t isa SymbolicIndex{Int}), keys(dict))
        if any(t -> !(t isa SymbolicIndex), keys(dict))
            throw(ArgumentError("property-dict must have VIndex/EIndex as keys!"))
        end
        resolved_dict = Dict{SymbolicIndex{Int}, Any}()
        for (k, v) in pairs(dict)
            resolvedidx = idxtype(k)(resolvecompidx(nw, k))
            resolved_dict[resolvedidx] = v
        end
        return resolved_dict
    else
        return dict
    end
end

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
