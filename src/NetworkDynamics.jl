module NetworkDynamics

include("ScalarVariables.jl") #not sure if we should keep that. Functionality can actually be easily incorporated into the other functions via kwargs.
using .ScalarVariables

include("StaticLines.jl")
using .StaticLines

include("DynamicLines.jl")
using .DynamicLines

end # module
