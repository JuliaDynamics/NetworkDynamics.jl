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

function initialization_problem(cf::T; t=NaN, verbose=true) where {T<:ComponentModel}
    hasinsym(cf) || throw(ArgumentError("Component model musst have `insym`!"))

    outfree_ms = Tuple((!).(map(s -> has_default(cf, s), sv)) for sv in outsym_normalized(cf))
    outfixs = Tuple(Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sv] for sv in outsym_normalized(cf))

    ufree_m = (!).(map(Base.Fix1(has_default, cf), sym(cf)))
    ufix = Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sym(cf)]

    infree_ms = Tuple((!).(map(s -> has_default(cf, s), sv)) for sv in insym_normalized(cf))
    infixs = Tuple(Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sv] for sv in insym_normalized(cf))

    pfree_m = (!).(map(Base.Fix1(has_default, cf), psym(cf)))
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
        throw(ArgumentError("Initialization problem underconstraint. $(Neqs) Equations and  $(Nfree) variables."))
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

    # check for missin guesses
    missing_guesses = Symbol[]
    uguess = map(freesym) do s
        if has_guess(cf, s)
            Float64(get_guess(cf, s))
        else
            push!(missing_guesses, s)
        end
    end

    for s in freesym
        if has_bounds(cf, s)
            @warn "Ignore bounds $(get_bounds(cf, s)) for $s. Not supported yet."
        end
    end

    isempty(missing_guesses) || throw(ArgumentError("Missing guesses for free variables $(missing_guesses)"))

    N = ForwardDiff.pickchunksize(Nfree)
    fz = let fg = compfg(cf),
        outcaches=map(d->DiffCache(zeros(d), N), outdim_normalized(cf)),
        ucache=DiffCache(zeros(dim(cf)), N),
        incaches=map(d->DiffCache(zeros(d), N), indim_normalized(cf)),
        pcache=DiffCache(zeros(pdim(cf))),
        t=t,
        outfree_ms=outfree_ms, ufree_m=ufree_m, infree_ms=infree_ms, pfree_m=pfree_m,
        outfixs=outfixs, ufix=ufix, infixs=infixs, pfix=pfix,
        unl_range_outs=unl_range_outs, unl_range_u=unl_range_u,
        unl_range_ins=unl_range_ins, unl_range_p=unl_range_p

        (dunl, unl, _) -> begin
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
                buf[mask] .= unl[range]
            end
            ubuf[ufree_m] .= unl[unl_range_u]
            for (buf, mask, range) in zip(inbufs, infree_ms, unl_range_ins)
                buf[mask] .= unl[range]
            end
            pbuf[pfree_m] .= unl[unl_range_p]

            # view into du buffer for the fg funtion
            @views dunl_fg = dunl[1:dim(cf)]
            # view into the output buffer for the outputs
            @views dunl_out = dunl[dim(cf)+1:end]

            # this fills the second half of the du buffer with the fixed and current outputs
            dunl_out .= RecursiveArrayTools.ArrayPartition(outbufs...)
            # execute fg to fill dunl and outputs
            fg(outbufs..., dunl, ubuf, inbufs..., pbuf, t)

            # calculate the residual for the second half ov the dunl buf, the outputs
            dunl_out .= dunl_out .- RecursiveArrayTools.ArrayPartition(outbufs...)
            nothing
        end
    end

    nlf = NonlinearFunction(fz; resid_prototype=zeros(Neqs), sys=SII.SymbolCache(freesym))

    if Neqs == Nfree
        verbose && @info "Initialization problem is fully constrained. Created NonlinearLeastSquaresProblem for $freesym"
    else
        verbose && @info "Initialization problem is overconstrained ($Nfree vars for $Neqs equations). Create NonlinearLeastSquaresProblem for $freesym."
    end
    NonlinearLeastSquaresProblem(nlf, uguess)
end

"""
    initialize_component!(cf::ComponentModel; verbose=true, kwargs...)

Initialize a `ComponentModel` by solving the corresponding `NonlinearLeastSquaresProblem`.
During initialization, everyting which has a `default` value (see [Metadata](@ref)) is considered
"fixed". All other variables are considered "free" and are solved for. The initial guess for each
variable depends on the `guess` value in the [Metadata](@ref).

The result is stored in the `ComponentModel` itself. The values of the free variables are stored
in the metadata field `init`.

The `kwargs` are passed to the nonlinear solver.
"""
function initialize_component!(cf; verbose=true, kwargs...)
    prob = initialization_problem(cf; verbose)

    if !isempty(prob.u0)
        sol = SciMLBase.solve(prob; kwargs...)

        if sol.prob isa NonlinearLeastSquaresProblem && sol.retcode == SciMLBase.ReturnCode.Stalled
            # https://github.com/SciML/NonlinearSolve.jl/issues/459
            res = LinearAlgebra.norm(sol.resid)
            @warn "Initialization for componend stalled with residual $(res)"
        elseif !SciMLBase.successful_retcode(sol.retcode)
            throw(ArgumentError("Initialization failed. Solver returned $(sol.retcode)"))
        end
        set_init!.(Ref(cf), SII.variable_symbols(sol), sol.u)
        resid = sol.resid
    else
        resid = init_residual(cf; recalc=true)
    end

    set_metadata!(cf, :init_residual, resid)

    verbose && @info "Initialization successful with residual $(LinearAlgebra.norm(sol.resid))"
    cf
end

function isinitialized(cf::ComponentModel)
    all(has_default_or_init(cf, s) for s in vcat(sym(cf), psym(cf)))
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
        p = Float64[get_default_or_init(cf, s) for s in psym(cf)]
        res = zeros(dim(cf) + sum(outdim_normalized(cf)))

        res_fg  = @views res[1:dim(cf)]
        res_out = @views res[dim(cf)+1:end]
        res_out .= RecursiveArrayTools.ArrayPartition(outs...)
        compfg(cf)(outs..., res_fg, u, ins..., p, t)
        res_out .= res_out .- RecursiveArrayTools.ArrayPartition(outs...)

        LinearAlgebra.norm(res)
    else
        LinearAlgebra.norm(get_metadata(cf, :init_residual))
    end
end
