module NetworkDynamics

include("nd_ODE_ODE_scalar.jl")
using .nd_ODE_ODE_scalar

include("nd_ODE_Static_scalar.jl")
using .nd_ODE_Static_scalar

include("nd_ODE_ODE.jl")
using .nd_ODE_ODE

include("nd_ODE_Static.jl")
using .nd_ODE_Static

end
