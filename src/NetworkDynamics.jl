module NetworkDynamics

include("nd_ODE_ODE_scalar.jl")
using .nd_ODE_ODE_scalar_mod

include("nd_ODE_Static_scalar.jl")
using .nd_ODE_Static_scalar_mod

include("nd_ODE_ODE.jl")
using .nd_ODE_ODE_mod

include("nd_ODE_Static.jl")
using .nd_ODE_Static_mod

end
