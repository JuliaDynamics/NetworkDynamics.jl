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

function initialization_problem(cf::T; t=NaN, apply_bound_transformation, verbose=true) where {T<:ComponentModel}
    hasinsym(cf) || throw(ArgumentError("Component model musst have `insym`!"))

    outfree_ms = Tuple((!).(map(s -> has_default(cf, s), sv)) for sv in outsym_normalized(cf))
    outfixs = Tuple(Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sv] for sv in outsym_normalized(cf))

    ufree_m = (!).(map(Base.Fix1(has_default, cf), sym(cf)))
    ufix = Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sym(cf)]

    infree_ms = Tuple((!).(map(s -> has_default(cf, s), sv)) for sv in insym_normalized(cf))
    infixs = Tuple(Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sv] for sv in insym_normalized(cf))

    is_free_p = s -> !is_unused(cf, s) && !has_default(cf, s) # free p are not unused and have no default
    pfree_m = map(is_free_p, psym(cf))
    pfix = Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in psym(cf)]

    # count free variables and equations
    Nfree = mapreduce(sum, +, outfree_ms) + sum(ufree_m) + mapreduce(sum, +, infree_ms) + sum(pfree_m)
    Neqs = dim(cf)  + mapreduce(length, +, outfree_ms)


    freesym = vcat(mapreduce((syms, map) -> syms[map], vcat, outsym_normalized(cf), outfree_ms),
                   sym(cf)[ufree_m],
                   mapreduce((syms, map) -> syms[map], vcat, insym_normalized(cf), outfree_ms),
                   psym(cf)[pfree_m])
    @assert length(freesym) == Nfree
    if Neqs < Nfree
        throw(ArgumentError("Initialization problem underconstraint. $(Neqs) Equations and  $(Nfree) variables: $freesym"))
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
    bounds = map(freesym) do sym
        if has_bounds(cf, sym)
            bound = get_bounds(cf, sym)
            if bound[1] >= 0 && bound[2] > bound[1]
                return :pos
            elseif bound[1] < bound[2]  && bound[2] <= 0
                return :neg
            end
        end
        :none
    end
    if !apply_bound_transformation || all(isequal(:none), bounds)
        boundT! = identity
        inv_boundT! = identity
    else
        if verbose
            idxs = findall(!isequal(:none), bounds)
            @info "Apply positivity/negativity conserving variable transformation on $(freesym[idxs]) to satisfy bounds."
        end
        boundT! = (u) -> begin
            for i in eachindex(u, bounds)
                if bounds[i] == :pos
                    u[i] = u[i]^2
                elseif bounds[i] == :neg
                    u[i] = -u[i]^2
                end
            end
            return u
        end
        inv_boundT! = (u) -> begin
            for i in eachindex(u, bounds)
                if bounds[i] == :pos
                    u[i] = sqrt(u[i])
                elseif bounds[i] == :neg
                    u[i] = sqrt(-u[i])
                end
            end
            return u
        end
    end

    # check for missin guesses
    missing_guesses = Symbol[]
    uguess = map(freesym) do s
        if has_guess(cf, s)
            Float64(get_guess(cf, s))
        else
            push!(missing_guesses, s)
        end
    end
    isempty(missing_guesses) || throw(ArgumentError("Missing guesses for free variables $(missing_guesses)"))

    # apply bound conserving transformation to initial state
    try
        inv_boundT!(uguess)
    catch
        throw(ArgumentError("Initial guess violates bounds. Check the docstring on `NetworkDynamics.initialize_component!`\
                             about bound satisfying transformations!"))
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
    initialize_component!(cf::ComponentModel; verbose=true, apply_bound_transformation=true, kwargs...)

Initialize a `ComponentModel` by solving the corresponding `NonlinearLeastSquaresProblem`.
During initialization, everyting which has a `default` value (see [Metadata](@ref)) is considered
"fixed". All other variables are considered "free" and are solved for. The initial guess for each
variable depends on the `guess` value in the [Metadata](@ref).

The result is stored in the `ComponentModel` itself. The values of the free variables are stored
in the metadata field `init`.

The `kwargs` are passed to the nonlinear solver.

## Bounds of free variables
When encountering any bounds in the free variables, NetworkDynamics will try to conserve them
by applying a coordinate transforamtion. This behavior can be supressed by setting `apply_bound_transformation`.
The following transformations are used:
- (a, b) intervals where both a and b are positive are transformed to `u^2`/`sqrt(u)`
- (a, b) intervals where both a and b are negative are transformed to `-u^2`/`sqrt(-u)`
"""
function initialize_component!(cf; verbose=true, apply_bound_transformation=true, kwargs...)
    prob, boundT! = initialization_problem(cf; verbose, apply_bound_transformation)

    if !isempty(prob.u0)
        sol = SciMLBase.solve(prob; verbose, kwargs...)

        if sol.prob isa NonlinearLeastSquaresProblem && sol.retcode == SciMLBase.ReturnCode.Stalled
            # https://github.com/SciML/NonlinearSolve.jl/issues/459
            res = LinearAlgebra.norm(sol.resid)
            @warn "Initialization for componend stalled with residual $(res)"
        elseif !SciMLBase.successful_retcode(sol.retcode)
            throw(ArgumentError("Initialization failed. Solver returned $(sol.retcode)"))
        else
            verbose && @info "Initialization successful with residual $(LinearAlgebra.norm(sol.resid))"
        end
        # transform back to original space
        u = boundT!(copy(sol.u))
        set_init!.(Ref(cf), SII.variable_symbols(sol), u)
        resid = sol.resid
    else
        resid = init_residual(cf; recalc=true)
        verbose && @info "No free variables! Residual $(LinearAlgebra.norm(resid))"
    end

    set_metadata!(cf, :init_residual, resid)

    broken = broken_bounds(cf)
    if !isempty(broken)
        @warn "Initialized model has broken bounds for $(broken). Use `dump_initial_state(mode)` \
        to inspect further and try to adapt the initial guesses!"
    end
    cf
end

function isinitialized(cf::ComponentModel)
    all(has_default_or_init(cf, s) || is_unused(cf, s) for s in vcat(sym(cf), psym(cf)))
end

"""
    init_residual(cf::T; t=NaN, recalc=false)

Calculates the residual |du| for the given component model for the values
provided via `default` and `init` [Metadata](@ref).

If recalc=false just return the residual determined in the actual initialization process.

See also [`initialize_component!`](@ref).
"""
function init_residual(cf::T; t=NaN, recalc=false) where {T<:ComponentModel}
    if !isinitialized(cf)
        throw(ArgumentError("Component is not initialized."))
    end

    if recalc || !has_metadata(cf, :init_residual)
        outs = Tuple(Float64[get_default_or_init(cf, s) for s in sv] for sv in outsym_normalized(cf))
        u = Float64[get_default_or_init(cf, s) for s in sym(cf)]
        ins = Tuple(Float64[get_default_or_init(cf, s) for s in sv] for sv in insym_normalized(cf))
        p = Float64[is_unused(cf, s) ? NaN : get_default_or_init(cf, s) for s in psym(cf)]
        res = zeros(dim(cf) + sum(outdim_normalized(cf)))

        res_fg  = @views res[1:dim(cf)]
        res_out = @views res[dim(cf)+1:end]
        res_out .= RecursiveArrayTools.ArrayPartition(outs...)
        compfg(cf)(outs, res_fg, u, ins, p, t)
        res_out .= res_out .- RecursiveArrayTools.ArrayPartition(outs...)

        LinearAlgebra.norm(res)
    else
        LinearAlgebra.norm(get_metadata(cf, :init_residual))
    end
end

function broken_bounds(cf)
    allsyms = vcat(sym(cf), psym(cf), insym_all(cf), outsym_flat(cf), obssym(cf))
    boundsyms = filter(s -> has_bounds(cf, s), allsyms)
    bounds = get_bounds.(Ref(cf), boundsyms)
    vals = get_initial_state(cf, boundsyms)
    broken = filter(i -> !bounds_satisfied(vals[i], bounds[i]), 1:length(bounds))
    boundsyms[broken]
end

function bounds_satisfied(val, bounds)
    @assert length(bounds) == 2
    !isnothing(val) && !isnan(val) && first(bounds) ≤ val ≤ last(bounds)
end

"""
    set_defaults!(nw::Network, s::NWState)

Set the default values of the network to the values of the given state.
Can be used to "store" the found fixpoint in the network metadata.
"""
function set_defaults!(nw::Network, s::NWState)
    for (sni, val) in zip(SII.variable_symbols(nw), uflat(s))
        set_default!(nw, sni, val)
    end
    for (sni, val) in zip(SII.parameter_symbols(nw), pflat(s))
        set_default!(nw, sni, val)
    end
    nw
end
