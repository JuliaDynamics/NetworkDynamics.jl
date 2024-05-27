using NDPrototype
using Graphs
using OrdinaryDiffEq
import SymbolicIndexingInterface as SII

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

g = complete_graph(3)
vf = Lib.kuramoto_second()
ef = Lib.kuramoto_edge()
nw = Network(g, vf, ef)
u0 = rand(dim(nw))
p = rand(pdim(nw))

odef = ODEFunction(nw; sys=nw)
prob = ODEProblem(odef, u0, (0,10), p)
sol = solve(prob, Tsit5())

using NDPrototype: VIndex
@test SII.is_variable(nw, VIndex(1,:δ))
@test !SII.is_variable(nw, VIndex(1,:F))
@test_throws BoundsError !SII.is_variable(nw, VIndex(10,:F))

@test sol[VIndex(1,:δ)] == sol[1,:]
@test sol[VIndex(1,:ω)] == sol[2,:]
@test sol[VIndex(3,:ω)] == sol[6,:]


@test !SII.is_variable(nw, VIndex(1:3,:ω))
@test collect(VIndex(1, [:δ,:ω])) == [VIndex(1,:δ), VIndex(1,:ω)]
@test collect(VIndex(1, [:δ,2])) == [VIndex(1,:δ), VIndex(1,2)]
@test collect(VIndex(1, 1:2)) == [VIndex(1,1), VIndex(1,2)]
@test collect(VIndex(1:5, 2)) == [VIndex(i,2) for i in 1:5]
@test collect(VIndex([1,3], :δ)) == [VIndex(1,:δ), VIndex(3,:δ)]
@test_throws ErrorException collect(VIndex(1:2, 1:5))
@test !SII.is_variable(nw, VIndex(3,[1,:ω]))

@test sol[VIndex(1:3,:ω),1] == u0[[2,4,6]]
@test sol[VIndex(3,[:δ,:ω]),1] == u0[[5,6]]
@test sol[VIndex(3,[:ω,:δ]),1] == u0[[6,5]]
@test sol[VIndex(3,[1,:ω]),1] == u0[[5,6]]

SII.variable_symbols(sol)
SII.parameter_symbols(sol)

SII.default_values(sol)
