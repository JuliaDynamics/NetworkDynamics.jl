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
#### more complex problem
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
Main.test_execution_styles(prob) # testing all ex styles #src
sol = solve(prob, Tsit5())

@test SII.variable_index.(Ref(nw), SII.variable_symbols(nw)) == 1:dim(nw)
@test SII.parameter_index.(Ref(nw), SII.parameter_symbols(nw)) == 1:pdim(nw)

####
#### State tests
####
using NetworkDynamics: NWState, NWParameter
t = 1.0
_uflat = copy(sol(t))
_pflat = copy(sol.prob.p)
s = NWState(nw, _uflat, _pflat)

SII.getu(s, EIndex(1,:e_dst))(s)
SII.getp(s, VPIndex(1,:M))(s)

SII.is_variable(nw, EIndex(1,:e_dst))
SII.variable_index(nw, EIndex(1,:e_dst))

@test map(idx->s[idx], SII.variable_symbols(nw)) == _uflat
@test map(idx->s[idx], SII.parameter_symbols(nw)) == _pflat
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

p = s.p
@test s.v[[1,3],:δ] == s[vidxs(s,:,"δ")]
@test_throws DimensionMismatch s.v[[1,3],:δ] = 1
@test_throws DimensionMismatch s.p.e[[1,5], :τ] = 1

s.p.e[[1,5], :τ] .= 1
@test s.p.e[[1,5], :τ] == [1,1]
s.v[[1,3],:δ] .= (1,2)
@test s.v[[1,3],:δ] == [1,2]
@test_throws DimensionMismatch s.v[1,:δ] = (1,2)
@test_throws DimensionMismatch s.p.e[1,:τ] = [1,2]
@test_throws DimensionMismatch s[vidxs(s,:,"δ")] = 1
@test s.p[vpidxs(s,:,"M")] == s.p[[VIndex(1,:M), VIndex(3,:M)]]

s.p[vpidxs(s,:,"M")] .= 4
@test s[vpidxs(s,:,"M")] == [4,4]
s[vpidxs(s,:,"M")] .= 5
@test s[vpidxs(s,:,"M")] == [5,5]

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
@test_broken s[EIndex(:,1)] == s[EIndex(1:6,1)]

          # variable     variable     observed     observed
mixedidx = [VIndex(1,1), EIndex(1,1), EIndex(2,1), EIndex(3,1)]
@test SII.is_observed.(s, mixedidx) == [0,0,1,1]
@test SII.is_variable.(s, mixedidx) == [1,1,0,0]
s[mixedidx] # -> calls observed for all indices

mixedidx = [VIndex(1,1), EIndex(1,1), EIndex(2,1)]
@test SII.is_observed.(s, mixedidx) == [0,0,1]
@test SII.is_variable.(s, mixedidx) == [1,1,0]
s[mixedidx] # -> calls observed for all indices

@test_broken s[EIndex(:,1)]
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
    # EIndex(1,:), # observed bc
    # EIndex(2,:), # variable bc
    VPIndex(1,:), # parameter bc
    EPIndex(1,:), # parameter bc
    VIndex(:,1), # variable bc first
    # EIndex(:,1), # observed bc first
    # EIndex(:,1), # variable bc first
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

@info "Test state getindex call"
for idx in idxtypes
    b = @b $s[$idx]
    if b.allocs != 0
        println(idx, " => ", b.allocs, " allocations")
    end
    @test b.allocs <= 16
end

# tests for state/parameter constructing/conversion
using NetworkDynamics: _init_flat, filltype
T = Vector{Float64}
@test isequal(_init_flat(T, 10, filltype(T)), [NaN for _ in 1:10])
T = Vector{Int64}
@test _init_flat(T, 10, filltype(T)) == zeros(10)
T = Vector{Union{Int64, Nothing}}
@test _init_flat(T, 10, filltype(T)) == [nothing for _ in 1:10]
T = Vector{Union{Int64, Missing}}
@test isequal(_init_flat(T, 10, filltype(T)), [missing for _ in 1:10])
@test isequal(pflat(NWState(nw)), pflat(NWParameter(nw)))
@test isequal(pflat(NWState(p; ufill=0)), pflat(p))

p = NWParameter(nw)
p.e[2:3,:K] .= 0
@test p.e[2:3,:K] == [0,0]
p.e[2:3,:K] = [0,0]
@test p.e[2:3,:K] == [0,0]

@test eltype(p) == Float64
p2 = NWParameter(p; ptype=Vector{Union{Float64,Nothing}})
@test eltype(p2) == Union{Float64,Nothing}

p3 = NWParameter(p)
p.v[1,:M] = 42
@test p3.v[1,:M] != 42


s = NWState(nw)
s[:] .= 1
s2 = NWState(s, utype=Vector{Float64})
@test eltype(s2) == Float64


s.p[:] .= 0
s3 = NWState(s, utype=Vector{Int}, ptype=Vector{Int})
@test eltype(NetworkDynamics.uflat(s3)) == Int
@test eltype(NetworkDynamics.uflat(s3)) == Int

# test new index generator methods
using NetworkDynamics: vidxs, eidxs, vpidxs, epidxs
n1 = ODEVertex(x->x, [:u, :v], [:p1, :p2])
n2 = ODEVertex(x->x, [:x1, :x2], [:p1, :p2], obsf=identity, obssym=[:o1, :o2])
n3 = ODEVertex(x->x, 3; name=:Vertex3)
e1 = StaticEdge(x->x, [:e1], AntiSymmetric())
e2 = StaticEdge(x->x, [:esrc,:edst], Fiducial())
g = path_graph(3)
nw = Network(g, [n1, n2, n3], [e1, e2])

@test vidxs(nw) == [VIndex(1, :u),
                    VIndex(1, :v),
                    VIndex(2, :x1),
                    VIndex(2, :x2),
                    VIndex(3, :v₁),
                    VIndex(3, :v₂),
                    VIndex(3, :v₃)]
@test eidxs(nw) == [EIndex(1, :e1),
                    EIndex(2, :esrc),
                    EIndex(2, :edst)]

@test vidxs(nw, 1) == [VIndex(1, :u), VIndex(1, :v)]
@test vidxs(nw, :ODEVertex) == [VIndex(1, :u),
                                VIndex(1, :v),
                                VIndex(2, :x1),
                                VIndex(2, :x2)]
@test vidxs(nw, "3") == [VIndex(3, :v₁),
                         VIndex(3, :v₂),
                         VIndex(3, :v₃)]
@test vidxs(nw, :, "v") == [VIndex(1, :v),
                            VIndex(3, :v₁),
                            VIndex(3, :v₂),
                            VIndex(3, :v₃)]
@test vidxs(nw, 2, "v") == VIndex[]
@test vpidxs(nw,:,"p") == vpidxs(nw,:,:) == [VPIndex(1, :p1),
                                             VPIndex(1, :p2),
                                             VPIndex(2, :p1),
                                             VPIndex(2, :p2)]
@test epidxs(nw,:,:) == EPIndex[]
