using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using Chairmarks
using Test
import SymbolicIndexingInterface as SII
using NetworkDynamics: VIndex, EIndex, VPIndex, EPIndex, _resolve_colon

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

@test SII.is_variable(nw, VIndex(1:3,:ω))
@test collect(VIndex(1, [:δ,:ω])) == [VIndex(1,:δ), VIndex(1,:ω)]
@test collect(VIndex(1, [:δ,2])) == [VIndex(1,:δ), VIndex(1,2)]
@test collect(VIndex(1, 1:2)) == [VIndex(1,1), VIndex(1,2)]
@test collect(VIndex(1:5, 2)) == [VIndex(i,2) for i in 1:5]
@test collect(VIndex([1,3], :δ)) == [VIndex(1,:δ), VIndex(3,:δ)]
@test_throws MethodError collect(VIndex(1:2, 1:5))
@test SII.is_variable(nw, VIndex(3,[1,:ω]))
@test !SII.is_variable(nw, VIndex(3,[1,:ω,:foobar]))

# test indexing of multipe variables
@test sol[VIndex(1:3,:ω),1] == u0[[2,4,6]]
@test sol[VIndex(3,[:δ,:ω]),1] == u0[[5,6]]
@test sol[VIndex(3,[:ω,:δ]),1] == u0[[6,5]]
@test sol[VIndex(3,[1,:ω]),1] == u0[[5,6]]

@test NetworkDynamics.observed_symbols(nw) == [EIndex(1, :P), EIndex(2, :P), EIndex(3, :P)]

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

@test filter(s->SII.is_observed(nw,s), SII.all_symbols(nw)) == NetworkDynamics.observed_symbols(nw)
@test filter(s->SII.is_parameter(nw,s), SII.all_symbols(nw)) == SII.parameter_symbols(nw)
@test filter(s->SII.is_variable(nw,s), SII.all_symbols(nw)) == SII.variable_symbols(nw)

sol[EIndex(1,:P)]
sol[EIndex(2,:P)]
sol[EIndex(3,:P)]
@test sol[EIndex(1:3,:P)] == sol[[EIndex(1,:P),EIndex(2,:P),EIndex(3,:P)]]

####
#### more complex prolbme
####
g = complete_graph(4)
vf = [Lib.kuramoto_second(), Lib.diffusion_vertex(), Lib.kuramoto_second(), Lib.diffusion_vertex()]
ef = [Lib.diffusion_odeedge(),
      Lib.kuramoto_edge(),
      Lib.kuramoto_edge(),
      Lib.diffusion_edge_fid(),
      Lib.diffusion_odeedge(),
      Lib.diffusion_edge_fid()]
nw = Network(g, vf, ef)
prob = ODEProblem(nw, rand(dim(nw)), (0,1), rand(pdim(nw)))
sol = solve(prob, Tsit5())

@test SII.variable_index.(Ref(nw), SII.variable_symbols(nw)) == 1:dim(nw)
@test SII.parameter_index.(Ref(nw), SII.parameter_symbols(nw)) == 1:pdim(nw)

####
#### State tests
####
using NetworkDynamics: NWState
t = 1.0
uflat = copy(sol(t))
pflat = copy(sol.prob.p)
s = NWState(nw, uflat, pflat)

SII.getu(s, EIndex(1,:e_dst))(s)
SII.getp(s, VPIndex(1,:M))(s)

SII.is_variable(nw, EIndex(1,:e_dst))
SII.variable_index(nw, EIndex(1,:e_dst))

@test map(idx->s[idx], SII.variable_symbols(nw)) == uflat
@test map(idx->s[idx], SII.parameter_symbols(nw)) == pflat
NetworkDynamics.observed_symbols(nw) .=> map(idx->s[idx], NetworkDynamics.observed_symbols(nw))

@test s[VPIndex(1,:M)] != 1.0
s[VPIndex(1,:M)] = 1
@test s[VPIndex(1,:M)] == 1.0

@test s[VIndex(1,2)] != 1.0
s[VIndex(1,2)] = 1.0
@test s[VIndex(1,2)] == 1.0

@test s.v[1,1] == s[VIndex(1,1)]
s.v[1,1] = 15
@test s.v[1,1] == s[VIndex(1,1)] == 15

@test s.e[5,:e_src] == s[EIndex(5,:e_src)]
s.v[1,1] = 15
@test s.v[1,1] == s[VIndex(1,1)] == 15

@test s.p.v[1,1] == s[VPIndex(1,1)]
s.p.v[1,1] = 10
@test s.p.v[1,1] == s[VPIndex(1,1)] == 10

@test s.p.e[1,1] == s[EPIndex(1,1)]
s.p.e[1,1] = 10
@test s.p.e[1,1] == s[EPIndex(1,1)] == 10

####
#### Tests for index with colon
####
@test NetworkDynamics._resolve_colon(nw, VIndex(:, :δ)) == VIndex(1:nv(g), :δ)
@test NetworkDynamics._resolve_colon(nw, EIndex(:, :δ)) == EIndex(1:ne(g), :δ)
@test NetworkDynamics._resolve_colon(nw, VPIndex(:, :δ)) == VPIndex(1:nv(g), :δ)
@test NetworkDynamics._resolve_colon(nw, EPIndex(:, :δ)) == EPIndex(1:ne(g), :δ)
@test NetworkDynamics._resolve_colon(nw, VIndex(3, :)) == VIndex(3, 1:2)
@test NetworkDynamics._resolve_colon(nw, VPIndex(3, :)) == VPIndex(3, 1:3)
@test NetworkDynamics._resolve_colon(nw, EIndex(4, :)) == EIndex(4, 1:2)
@test NetworkDynamics._resolve_colon(nw, EPIndex(4, :)) == EPIndex(4, 1:1)

@inferred NetworkDynamics._resolve_colon(nw, VIndex(:, :δ))
@inferred NetworkDynamics._resolve_colon(nw, EIndex(:, :δ))
@inferred NetworkDynamics._resolve_colon(nw, VPIndex(:, :δ))
@inferred NetworkDynamics._resolve_colon(nw, EPIndex(:, :δ))
@inferred NetworkDynamics._resolve_colon(nw, VIndex(3, :))
@inferred NetworkDynamics._resolve_colon(nw, VPIndex(3, :))
@inferred NetworkDynamics._resolve_colon(nw, EIndex(4, :))
@inferred NetworkDynamics._resolve_colon(nw, EPIndex(4, :))

for et in [VIndex, EIndex, VPIndex, EPIndex]
    collect(et(1:2,:δ))
    collect(et(1,1:5))
    collect(et(1,[:foo, :bar]))
    repr.(et(1:2,:δ))
    repr.(et(1,:δ))
end

@test s[[VIndex(1,1), VPIndex(1,2)]] == [s[VIndex(1,1)], s[VPIndex(1,2)]]
@test s[(VIndex(1,1), VPIndex(1,2))] == (s[VIndex(1,1)], s[VPIndex(1,2)])
@test s[[VIndex(1,1), VIndex(1,2)]] == [s[VIndex(1,1)], s[VIndex(1,2)]]
@test s[(VIndex(1,1), VIndex(1,2))] == (s[VIndex(1,1)], s[VIndex(1,2)])

@test s[VIndex(:,1)] == s[VIndex(1:4,1)]
@test s[EIndex(:,1)] == s[EIndex(1:6,1)]

          # variable     variable     observed     observed
mixedidx = [VIndex(1,1), EIndex(1,1), EIndex(2,1), EIndex(3,1)]
@test SII.is_observed.(s, mixedidx) == [0,0,1,1]
@test SII.is_variable.(s, mixedidx) == [1,1,0,0]
s[mixedidx] # -> calls observed for all indices

mixedidx = [VIndex(1,1), EIndex(1,1), EIndex(2,1)]
@test SII.is_observed.(s, mixedidx) == [0,0,1]
@test SII.is_variable.(s, mixedidx) == [1,1,0]
s[mixedidx] # -> calls observed for all indices

s[EIndex(:,1)]
idx = EIndex(:,1)

@test !SII.is_variable(s, EIndex(:,1))

using NetworkDynamics: _expand_and_collect
@test _expand_and_collect(nw, VIndex(1,1)) == VIndex(1,1)
@test _expand_and_collect(nw, VIndex(:,1)) == collect(VIndex(1:nv(g),1))
@test _expand_and_collect(nw, [VIndex(:,1)]) == collect(VIndex(1:nv(g),1))
@test _expand_and_collect(nw, [EIndex(4,5), VIndex(:,1)]) == vcat(EIndex(4,5), collect(VIndex(1:nv(g),1)))
@test _expand_and_collect(nw, [EIndex(1:2,5), VIndex(:,1)]) == vcat(collect(EIndex(1:2,5)), collect(VIndex(1:nv(g),1)))

####
#### check on performance
####

idxtypes = [
    VIndex(1,1), # variable
    EIndex(1,1), # observed
    EIndex(2,1), # variable
    VPIndex(1,1), # parameter
    EPIndex(1,1), # parameter
    VIndex(1,1:1), # variable bc
    EIndex(1,1:1), # observed bc
    EIndex(2,1:1), # variable bc
    VPIndex(1,1:1), # parameter bc
    EPIndex(1,1:1), # parameter bc
    VIndex(1,:), # variable bc
    EIndex(1,:), # observed bc
    EIndex(2,:), # variable bc
    VPIndex(1,:), # parameter bc
    EPIndex(1,:), # parameter bc
    VIndex(:,1), # variable bc first
    EIndex(:,1), # observed bc first
    EIndex(:,1), # variable bc first
    # VPIndex(:,1), # parameter bc first
    EPIndex(:,1), # parameter bc first
    VIndex(1:3,1), # variable bc first
    EIndex(1:3,1), # observed bc first
    EIndex(1:3,1), # variable bc first
    VPIndex([1,3],1), # parameter bc first
    EPIndex(1:3,1), # parameter bc first
]

using NetworkDynamics: _is_variable
for idx in idxtypes
    println("Test $idx")
    s[idx]
    # @inferred SII.is_variable(s, idx)
end

# is_xxx methods
@info "Test is_variable"
for idx in idxtypes
    b = @b SII.is_variable($s,$idx)
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs==0
end
@info "Test is_parameter"
for idx in idxtypes
    b = @b SII.is_parameter($s,$idx)
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs==0
end
@info "Test is_observed"
for idx in idxtypes
    b = @b SII.is_observed($s,$idx)
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs==0
end

# index methods
@info "Test variable_index"
for idx in idxtypes
    SII.is_variable(s, idx) || continue
    b = @b SII.variable_index($s,$idx)
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs <= 2 # 2 are used to create an array
end
@info "Test parameter_index"
for idx in idxtypes
    SII.is_parameter(s, idx) || continue
    b = @b SII.parameter_index($s,$idx)
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs <= 2 # 2 are used to create an array
end
@info "Test observed"
for idx in idxtypes
    SII.is_observed(s, idx) || continue
    b = @b SII.observed($s,$idx)
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs <= 4
end

@info "Test full call"
for idx in idxtypes
    b = @b $s[$idx]
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs <= 11
end
