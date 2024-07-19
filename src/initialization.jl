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
    _prob = NonlinearProblem(nw, prob)
    SciMLBase.solve(_prob, alg; kwargs...)
end

function _solve_fixpoint(prob, alg::SteadyStateDiffEqAlgorithm; kwargs...)
    sol = SciMLBase.solve(prob, alg; kwargs...)
end
