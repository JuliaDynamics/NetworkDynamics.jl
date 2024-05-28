using NDPrototype
using Graphs
using OrdinaryDiffEq
using Test
import SymbolicIndexingInterface as SII
using NDPrototype: VIndex, EIndex, VPIndex, EPIndex

(isinteractive() && @__MODULE__()==Main ? includet : include)("ComponentLibrary.jl")

g = complete_graph(3)
vf = Lib.kuramoto_second()
ef = Lib.kuramoto_edge()
nw = Network(g, vf, ef)
u0 = rand(dim(nw))
p = rand(pdim(nw))

SII.default_values(nw)
prob = ODEProblem(nw, u0, (0,10), p)
sol = solve(prob, Tsit5())

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
SII.default_values(sol)

@test NDPrototype.observed_symbols(nw) == [EIndex(1, :P), EIndex(2, :P), EIndex(3, :P)]

@test SII.parameter_symbols(nw) == [VPIndex{Int64, Symbol}(1, :M),
                                    VPIndex{Int64, Symbol}(1, :D),
                                    VPIndex{Int64, Symbol}(1, :Pm),
                                    VPIndex{Int64, Symbol}(2, :M),
                                    VPIndex{Int64, Symbol}(2, :D),
                                    VPIndex{Int64, Symbol}(2, :Pm),
                                    VPIndex{Int64, Symbol}(3, :M),
                                    VPIndex{Int64, Symbol}(3, :D),
                                    VPIndex{Int64, Symbol}(3, :Pm),
                                    EPIndex{Int64, Symbol}(1, :K),
                                    EPIndex{Int64, Symbol}(2, :K),
                                    EPIndex{Int64, Symbol}(3, :K)]
SII.all_variable_symbols(nw)

@test filter(s->SII.is_observed(nw,s), SII.all_symbols(nw)) == NDPrototype.observed_symbols(nw)
@test filter(s->SII.is_parameter(nw,s), SII.all_symbols(nw)) == SII.parameter_symbols(nw)
@test filter(s->SII.is_variable(nw,s), SII.all_symbols(nw)) == SII.variable_symbols(nw)

sol[EIndex(1,:P)]
sol[EIndex(2,:P)]
sol[EIndex(3,:P)]
@test sol[EIndex(1:3,:P)] == sol[[EIndex(1,:P),EIndex(2,:P),EIndex(3,:P)]]

g = complete_graph(4)
vf = [Lib.kuramoto_second(), Lib.diffusion_vertex(), Lib.kuramoto_second(), Lib.diffusion_vertex()]
ef = [Lib.diffusion_odeedge(),
      Lib.kuramoto_edge(),
      Lib.kuramoto_edge(),
      Lib.diffusion_edge_fid(),
      Lib.diffusion_odeedge(),
      Lib.diffusion_edge_fid()]
nw = Network(g, vf, ef)

@test SII.variable_index.(Ref(nw), SII.variable_symbols(nw)) == 1:dim(nw)
@test SII.parameter_index.(Ref(nw), SII.parameter_symbols(nw)) == 1:pdim(nw)
