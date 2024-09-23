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

function initialization_problem(cf::T; t=NaN, verbose=true) where {T<:Union{ODEVertex, ODEEdge}}
    ufree_m = (!).(map(Base.Fix1(has_default, cf), sym(cf)))
    pfree_m = (!).(map(Base.Fix1(has_default, cf), psym(cf)))
    ufree = sum(ufree_m)
    pfree = sum(pfree_m)
    ufix = Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in sym(cf)]
    pfix = Float64[has_default(cf, s) ? get_default(cf, s) : NaN for s in psym(cf)]

    hasinputsym(cf) || throw(ArgumentError("Vertex musst have `inputsym` with default values!"))
    input = if T <: EdgeFunction
        (;src=Float64[get_default(cf, s) for s in inputsym(cf).src], dst=Float64[get_default(cf, s) for s in inputsym(cf).dst])
    else
        Float64[get_default(cf, s) for s in inputsym(cf)]
    end

    freesym = vcat(sym(cf)[ufree_m], psym(cf)[pfree_m])

    Nfree = length(freesym)
    Neqs = dim(cf)
    if Neqs < Nfree
        throw(ArgumentError("Initialization problem underconstraint. $(Neqs) Equations and  $(Nfree) variables."))
    end

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

    fz = let f = compf(cf),
        ucache=DiffCache(zeros(dim(cf))),
        pcache=DiffCache(zeros(pdim(cf))),
        t=t,
        ufree_m=ufree_m, pfree_m=pfree_m, ufree=ufree, pfree=pfree,
        ufix=ufix, pfix=pfix,
        input=input,
        T=T

        (du, u, _) -> begin
            ubuf = PreallocationTools.get_tmp(ucache, u)
            pbuf = PreallocationTools.get_tmp(pcache, u)

            ubuf .= ufix
            ubuf[ufree_m] .= u[1:ufree]
            pbuf .= pfix
            pbuf[pfree_m] .= u[ufree+1:ufree+pfree]

            if T <: ODEEdge
                f(du, ubuf, input.src, input.dst, pbuf, t)
            elseif T <: ODEVertex
                f(du, ubuf, input, pbuf, t)
            end
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
    initialize_component!(cf::ComponentFunction; verbose=true, kwargs...)

Initialize a `ComponentFunction` by solving the corresponding `NonlinearLeastSquaresProblem`.
During initialization, everyting which has a `default` value (see [Metadata](@ref)) is considered
"fixed". All other variables are considered "free" and are solved for. The initial guess for each
variable depends on the `guess` value in the [Metadata](@ref).

The result is stored in the `ComponentFunction` itself. The values of the free variables are stored
in the metadata field `init`.

The `kwargs` are passed to the nonlinear solver.
"""
function initialize_component!(cf; verbose=true, kwargs...)
    prob = initialization_problem(cf)
    sol = SciMLBase.solve(prob; kwargs...)

    if sol.prob isa NonlinearLeastSquaresProblem && sol.retcode == SciMLBase.ReturnCode.Stalled
        # https://github.com/SciML/NonlinearSolve.jl/issues/459
        res = LinearAlgebra.norm(sol.resid)
        @warn "Initialization for componend stalled with residual $(res)"
    elseif !SciMLBase.successful_retcode(sol.retcode)
        throw(ArgumentError("Initialization failed. Solver returned $(sol.retcode)"))
    end
    set_init!.(Ref(cf), SII.variable_symbols(sol), sol.u)

    set_metadata!(cf, :init_residual, sol.resid)

    verbose && @info "Initialization successful with residual $(LinearAlgebra.norm(sol.resid))"
    cf
end

isinitialized(cf::StaticEdge) = true
isinitialized(cf::StaticVertex) = true
function isinitialized(cf::ComponentFunction)
    all(has_default_or_init(cf, s) for s in vcat(sym(cf), psym(cf)))
end

"""
    init_residual(cf::T; t=NaN, recalc=false)

Calculates the residual |du| for the given component function for the values
provided via `default` and `init` [Metadata](@ref).

If recalc=false just return the residual determined in the actual initialization process.

See also [`initialize_component!`](@ref).
"""
function init_residual(cf::T; t=NaN, recalc=false) where {T<:Union{ODEVertex, ODEEdge}}
    if !isinitialized(cf)
        throw(ArgumentError("Component is not initialized."))
    end

    if recalc || !has_metadata(cf, :init_residual)
        u = Float64[get_default_or_init(cf, s) for s in sym(cf)]
        p = Float64[get_default_or_init(cf, s) for s in psym(cf)]
        res = zeros(dim(cf))

        if T <: EdgeFunction
            src=Float64[get_default(cf, s) for s in inputsym(cf).src]
            dst=Float64[get_default(cf, s) for s in inputsym(cf).dst]
            compf(cf)(res, u, src, dst, p, t)
        else
            input = Float64[get_default(cf, s) for s in inputsym(cf)]
            compf(cf)(res, u, input, p, t)
        end
        LinearAlgebra.norm(res)
    else
        LinearAlgebra.norm(get_metadata(cf, :init_residual))
    end
end
