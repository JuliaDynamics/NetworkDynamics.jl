module Utilities

using NLsolve
using DifferentialEquations
using LinearAlgebra

export RootRhs
export find_valid_ic

struct RootRhs
    rhs
    mpm
end
function (rr::RootRhs)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    rr.mpm * dx .- dx
end

function RootRhs(of::ODEFunction)
    mm = of.mass_matrix
    @assert mm != nothing
    mpm = pinv(mm) * mm
    RootRhs(of.f, mpm)
end

function find_valid_ic(of::ODEFunction, ic_guess)
    rr = RootRhs(of)
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return nl_res.zero
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end

end #module
