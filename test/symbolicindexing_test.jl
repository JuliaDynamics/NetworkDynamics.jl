using NetworkDynamics
using Graphs
using OrdinaryDiffEqTsit5
using Chairmarks
using Test
using Symbolics
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

@test NetworkDynamics.observed_symbols(nw) == [EIndex(1, :₋P), EIndex(1, :P), EIndex(2, :₋P), EIndex(2, :P), EIndex(3, :₋P), EIndex(3, :P)]

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

# sol[obs] does not work, because obs has two timeseries: Continuous and Discrete
# it is unclear, whether it should return onlye discre values or for all sol.t
@test_broken sol[EIndex(1,:P)]
@test_broken sol[EIndex(2,:P)]
@test_broken sol[EIndex(3,:P)]
@test_broken sol[EIndex(1:3,:P)] == sol[[EIndex(1,:P),EIndex(2,:P),EIndex(3,:P)]]
@test sol(sol.t, idxs=EIndex(1:3,:P)) == sol(sol.t, idxs=[EIndex(1,:P),EIndex(2,:P),EIndex(3,:P)])
@test sol(sol.t, idxs=EIndex(1:3,:₋P)) == sol(sol.t, idxs=[EIndex(1,:₋P),EIndex(2,:₋P),EIndex(3,:₋P)])
@test sol(sol.t, idxs=EIndex(1,:P)).u == -1* sol(sol.t, idxs=EIndex(1,:₋P)).u


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

@test SII.getu(s, EIndex(1,:e_dst))(s) == uflat(s)[7]
@test SII.getp(s, VPIndex(1,:M))(s) ==pflat(s)[1]
@test SII.getp(s, VIndex(1,:M))(s) ==pflat(s)[1]

@test SII.is_variable(nw, EIndex(1,:e_dst))
@test SII.variable_index(nw, EIndex(1,:e_dst)) == 7

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
@test NetworkDynamics._resolve_colon(nw, EIndex(4, :)) == EIndex(4, 1:0)
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

@test s[[VIndex(1,1), VPIndex(1,2)]] == [s[VIndex(1,1)], s[VPIndex(1,2)]] == [uflat(s)[1], pflat(s)[2]]
@test s[(VIndex(1,1), VPIndex(1,2))] == (s[VIndex(1,1)], s[VPIndex(1,2)]) == (uflat(s)[1], pflat(s)[2])
@test s[[VIndex(1,1), VIndex(1,2)]] == [s[VIndex(1,1)], s[VIndex(1,2)]] == [uflat(s)[1], uflat(s)[2]]
@test s[(VIndex(1,1), VIndex(1,2))] == (s[VIndex(1,1)], s[VIndex(1,2)]) == (uflat(s)[1], uflat(s)[2])

@test s[VIndex(:,1)] == s[VIndex(1:4,1)]
@test_broken s[EIndex(:,1)] == s[EIndex(1:6,1)]

          # variable     variable     observed     observed
mixedidx = [VIndex(1,1), EIndex(1,1), EIndex(2,:P), EIndex(3,:P)]
@test SII.is_observed.(s, mixedidx) == [0,0,1,1]
@test SII.is_variable.(s, mixedidx) == [1,1,0,0]
s[mixedidx] # -> calls observed for all indices

mixedidx = [VIndex(1,1), EIndex(1,1), EIndex(2,:P)]
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
    EIndex(1,1), # variable
    EIndex(2,:P), # observed
    VPIndex(1,1), # parameter
    EPIndex(1,1), # parameter
    VIndex(1,1:1), # variable bc
    EIndex(1,1:1), # variable bc
    EIndex(2,[:P,:₋P]), # observed bc # does not work with ranges anymore
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
    EIndex(2:3,:P), # variable bc first # observables require symbols now
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
    if VERSION ≥ v"1.11"
        @test b.allocs <= 7
    end
end

@info "Test state getindex call"
for idx in idxtypes
    # HACK: _expand_and_collect to workaround https://github.com/SciML/SymbolicIndexingInterface.jl/issues/94
    _idx = NetworkDynamics._expand_and_collect(s, idx)
    getter = SII.getu(s, _idx)
    b = @b $(SII.getu)($s, $_idx)
    if b.allocs != 0
        println(rpad(idx,21), "=> ", b.allocs, " allocations to generate getter")
    end
    b = @b $getter($s)
    v = getter(s)
    if v isa Number
        @test b.allocs <= 0
        b.allocs != 0 && println(idx, " => ", b.allocs, " allocations to call getter")
    elseif v isa AbstractArray
        @test b.allocs <= 2
        b.allocs > 2 && println(idx, " => ", b.allocs, " allocations to call getter")
    else
        @test false
    end
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
fv = (du, u, ein, p, t) -> nothing
n1 = VertexModel(; f=fv, g=1:2, sym=[:u, :v], psym=[:p1, :p2], name=:VF)
n2 = VertexModel(; f=fv, g=1:2, sym=[:x1, :x2], psym=[:p1, :p2], obsf=identity, obssym=[:obs1, :obs2], name=:VF)
n3 = VertexModel(; f=fv, g=1:2, dim=3, name=:Vertex3)

gss = (odst, vsrc, vdst, p, t) -> nothing
gfid = (osrc, odst, vsrc, vdst, p ,t) -> nothing
e1 = EdgeModel(; g=AntiSymmetric(gss), outsym=[:e1])
e2 = EdgeModel(; g=gfid, outsym=(src=:esrc, dst=:edst))
g = path_graph(3)
nw = Network(g, [n1, n2, n3], [e1, e2])

@test vidxs(nw) == [VIndex(1, :u),
                    VIndex(1, :v),
                    VIndex(2, :x1),
                    VIndex(2, :x2),
                    VIndex(2, :obs1),
                    VIndex(2, :obs2),
                    VIndex(3, :v₁),
                    VIndex(3, :v₂),
                    VIndex(3, :v₃)]
@test eidxs(nw) == [EIndex(1, :₋e1),
                    EIndex(1, :e1),
                    EIndex(2, :esrc),
                    EIndex(2, :edst)]

@test vidxs(nw, 1) == [VIndex(1, :u), VIndex(1, :v)]
@test vidxs(nw, :VF) == [VIndex(1, :u),
                         VIndex(1, :v),
                         VIndex(2, :x1),
                         VIndex(2, :x2),
                         VIndex(2, :obs1),
                         VIndex(2, :obs2)]
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

####
#### Timeseries parameter test
####
using NetworkDynamics: vidxs, eidxs, vpidxs, epidxs
fv = (du, u, ein, p, t) -> nothing
n1 = VertexModel(; f=fv, g=1:2, sym=[:u, :v], psym=[:p1, :p2])
n2 = VertexModel(; f=fv, g=1:2, sym=[:x1, :x2], psym=[:p1, :p2], obsf=identity, obssym=[:obs1, :obs2])
n3 = VertexModel(; f=fv, g=1:2, dim=3, name=:Vertex3)

gss = (odst, vsrc, vdst, p, t) -> nothing
gfid = (osrc, odst, vsrc, vdst, p ,t) -> nothing
e1 = EdgeModel(; g=AntiSymmetric(gss), outsym=[:e1])
e2 = EdgeModel(; g=gfid, outsym=(src=:esrc, dst=:edst))
g = path_graph(3)
nw = Network(g, [n1, n2, n3], [e1, e2])

@test SII.get_all_timeseries_indexes(nw, :t) == Set([SII.ContinuousTimeseries()])
# https://github.com/SciML/SymbolicIndexingInterface.jl/issues/95
@test_broken SII.get_all_timeseries_indexes(nw, :x) == Set()

@test SII.get_all_timeseries_indexes(nw, VIndex(1,:u)) == Set([SII.ContinuousTimeseries()])
@test SII.get_all_timeseries_indexes(nw, VPIndex(1,:p1)) == Set([1])
@test SII.get_all_timeseries_indexes(nw, [VIndex(1,:u), VPIndex(1,:p1)]) == Set([SII.ContinuousTimeseries(), 1])

# test named vertices and edges
@testset "test sym indices for named edges/vertices" begin
    v1 = Lib.kuramoto_second(name=:v1)
    v2 = Lib.kuramoto_second(name=:v2)
    v3 = Lib.kuramoto_second(name=:v3)
    e1 = Lib.kuramoto_edge(name=:e1)
    e2 = Lib.kuramoto_edge(name=:e2)
    e3 = Lib.kuramoto_edge(name=:e3)
    g = complete_graph(3)
    nw = Network(g, [v1, v2, v3], [e1, e2, e3])
    s = NWState(nw, collect(1:dim(nw)), collect(dim(nw)+1:dim(nw)+pdim(nw)))
    @test_throws ArgumentError s.v[:v, 1]
    @test s.v[:v1, 1] == s[VIndex(1,1)]
    @test s.v[:v2, 1] == s[VIndex(2,1)]
    @test s.v[:v3, 1] == s[VIndex(3,1)]
    @test s.e[:e1, :P] == s[EIndex(1,:P)]
    @test s.e[:e2, :P] == s[EIndex(2,:P)]
    @test s.e[:e3, :P] == s[EIndex(3,:P)]
    @test_throws ArgumentError s.p.v[:v, 1]
    @test s.p.v[:v1, 1] == s[VPIndex(1,1)]
    @test s.p.v[:v2, 1] == s[VPIndex(2,1)]
    @test s.p.v[:v3, 1] == s[VPIndex(3,1)]
    @test s.p.e[:e1, 1] == s[EPIndex(1,1)]
    @test s.p.e[:e2, 1] == s[EPIndex(2,1)]
    @test s.p.e[:e3, 1] == s[EPIndex(3,1)]
end

# test observed for inputs
@testset "test observing of model input" begin
    v1 = Lib.kuramoto_second(name=:v1, vidx=1, insym=[:Pin])
    v2 = Lib.kuramoto_second(name=:v2, vidx=2, insym=[:Pin])
    v3 = Lib.kuramoto_second(name=:v3, vidx=3, insym=[:Pin])
    e1 = Lib.kuramoto_edge(name=:e1, src=1, dst=2, insym=[:δin])
    e2 = Lib.kuramoto_edge(name=:e2, src=2, dst=3, insym=[:δin])
    nw = Network([v1,v2,v3], [e1,e2])
    s = NWState(nw, rand(dim(nw)), rand(pdim(nw)))
    @test s[VIndex(:v1, :Pin)] == s[EIndex(:e1, :₋P)]
    @test s[VIndex(:v2, :Pin)] == s[EIndex(:e1, :P)] + s[EIndex(:e2, :₋P)]
    @test s[VIndex(:v3, :Pin)] == s[EIndex(:e2, :P)]
    @test s[EIndex(:e1, :src₊δin)] == s[VIndex(:v1, :δ)]
    @test s[EIndex(:e1, :dst₊δin)] == s[VIndex(:v2, :δ)]
    @test s[EIndex(:e2, :src₊δin)] == s[VIndex(:v2, :δ)]
    @test s[EIndex(:e2, :dst₊δin)] == s[VIndex(:v3, :δ)]
end

@testset "test observed expressions" begin
    v1 = Lib.kuramoto_second(name=:v1, vidx=1, insym=[:Pin])
    v2 = Lib.kuramoto_second(name=:v2, vidx=2, insym=[:Pin])
    v3 = Lib.kuramoto_second(name=:v3, vidx=3, insym=[:Pin])
    e1 = Lib.kuramoto_edge(name=:e1, src=1, dst=2, insym=[:δin])
    e2 = Lib.kuramoto_edge(name=:e2, src=2, dst=3, insym=[:δin])
    nw = Network([v1,v2,v3], [e1,e2])
    s = NWState(nw, rand(dim(nw)), rand(pdim(nw)))

    obsex = @obsex(VIndex(1,:δ) + VIndex(2,:δ))
    @test s[obsex] == s[VIndex(1,:δ)] + s[VIndex(2,:δ)]
    @test s[@obsex VIndex(1,:δ) - EIndex(:e1, :src₊δin)] == 0

    @test SII.getname(@obsex(VIndex(1,:δ) + VIndex(2,:δ))) == Symbol("v1₊δ+v2₊δ")
    @test SII.getname(@obsex(δ²=VIndex(1,:δ)^2)) == :δ²

    @test s[@obsex(vidxs(s, :, :δ) .- VIndex(1, :δ))] == [0, s[VIndex(2,:δ)] - s[VIndex(1,:δ)], s[VIndex(3,:δ)] - s[VIndex(1,:δ)]]

    obsex = @obsex(δ_rel = vidxs(s, :, :δ) .- VIndex(1, :δ))
end

@testset "test performace of created observed functions" begin
    v1 = Lib.kuramoto_second(name=:v1, vidx=1, insym=[:Pin])
    v2 = Lib.swing_mtk(name=:v2, vidx=2)
    set_default!(v2, :Pmech, -1.0)
    e = Lib.kuramoto_edge(name=:e12, src=1, dst=2)
    set_default!(e, :K, 1.0)
    nw = Network([v1,v2],e)

    s = NWState(nw)
    # normal state, observed and output state
    idxs1 = [VIndex(1,:δ), VIndex(2, :Pdamping), EIndex(1,:P), VIndex(2,:P)]
    idxs2 = [VIndex(1,:δ), VIndex(2,:θ)]
    # full call
    # @b $s[$idxs1] # 134 106 94 101
    # @b $s[$idxs2] # 31  31 34

    # scalar call
    # @b $s[$(VIndex(2,:Pdamping))] # 28

    # @b SII.observed($nw, $(VIndex(2,:Pdamping))) # 15
    # @b SII.observed($nw, $(VIndex(2,:θ))) # 7 5

    b = @b SII.observed($nw, $idxs1) # 69 36 42 30 37 47
    if VERSION ≥ v"1.11"
        @test b.allocs <= 47
    end
    b = @b SII.observed($nw, $idxs2) # 12 7 10 5
    if VERSION ≥ v"1.11"
        @test b.allocs <= 5
    end

    obsf1 = SII.observed(nw, idxs1)
    obsf2 = SII.observed(nw, idxs2)
    # @b $obsf1($(rand(dim(nw))), $(rand(pdim(nw))), NaN) # 81ns 2 allocs
    # @b $obsf2($(rand(dim(nw))), $(rand(pdim(nw))), NaN) # 34ns 2 allocs
    b = @b $obsf1($(rand(dim(nw))), $(rand(pdim(nw))), NaN, $(zeros(length(idxs1)))) # 64ns 0 allocs
    @test b.allocs == 0
    b = @b $obsf2($(rand(dim(nw))), $(rand(pdim(nw))), NaN, $(zeros(length(idxs2)))) # 17ns 0 allocs
    @test b.allocs == 0
end

@testset "test edge indexing with Pair syntax" begin
    # Create a simple network with named vertices for testing
    v1 = Lib.kuramoto_second(name=:v1)
    v2 = Lib.kuramoto_second(name=:v2)
    v3 = Lib.kuramoto_second(name=:v3)
    e12 = Lib.kuramoto_edge(name=:e12)
    e23 = Lib.kuramoto_edge(name=:e23)
    g = path_graph(3)
    nw = Network(g, [v1, v2, v3], [e12, e23])

    # Basic pair syntax tests
    @test SII.is_observed(nw, EIndex(1=>2, :P))
    @test SII.is_observed(nw, EIndex(:v1=>:v2, :P))
    @test SII.is_parameter(nw, EPIndex(1=>2, :K))
    @test SII.is_parameter(nw, EPIndex(:v1=>:v2, :K))

    # Test that pairs resolve to regular indices
    @test SII.parameter_index(nw, EPIndex(1=>2, :K)) == SII.parameter_index(nw, EPIndex(1, :K))
    @test SII.parameter_index(nw, EPIndex(:v1=>:v2, :K)) == SII.parameter_index(nw, EPIndex(1, :K))
    @test SII.parameter_index(nw, EPIndex(:v1=>2, :K)) == SII.parameter_index(nw, EPIndex(1, :K))

    # Error cases
    @test_throws ArgumentError SII.parameter_index(nw, EPIndex(1=>3, :K))  # no direct edge
    @test_throws ArgumentError SII.parameter_index(nw, EPIndex(:nonexistent=>:v2, :K))  # invalid vertex

    # Test with solution and state objects
    u0 = rand(dim(nw))
    p = rand(pdim(nw))
    prob = ODEProblem(nw, u0, (0.0, 1.0), p)
    sol = solve(prob, Tsit5())
    s = NWState(nw, u0, p)

    # Solution indexing
    @test sol([0.1, 0.5], idxs=EIndex(1=>2, :P)).u ≈ sol([0.1, 0.5], idxs=EIndex(1, :P)).u
    @test sol([0.1, 0.5], idxs=EIndex(:v1=>:v2, :P)).u ≈ sol([0.1, 0.5], idxs=EIndex(1, :P)).u

    # State access
    @test s[EIndex(1=>2, :P)] == s[EIndex(1, :P)]
    @test s[EIndex(:v1=>:v2, :P)] == s[EIndex(1, :P)]
    @test s[EPIndex(1=>2, :K)] == s[EPIndex(1, :K)]
    @test s[EPIndex(:v1=>:v2, :K)] == s[EPIndex(1, :K)]

    # Proxy syntax
    @test s.e[1=>2, :P] == s[EIndex(1, :P)]
    @test s.e[:v1=>:v2, :P] == s[EIndex(1, :P)]
    s.p.e[1=>2, :K] = 2.71
    @test s.p.e[1, :K] == 2.71
    @test s.p[EPIndex(1=>2, :K)] == 2.71

    # Multiple edge access
    edges_int = [EIndex(1=>2, :P), EIndex(2=>3, :P)]
    edges_named = [EIndex(:v1=>:v2, :P), EIndex(:v2=>:v3, :P)]
    edges_regular = [EIndex(1, :P), EIndex(2, :P)]
    @test s[edges_int] == s[edges_regular]
    @test s[edges_named] == s[edges_regular]

    # naming
    @test SII.getname(EIndex(1=>2, :P)) == :v1ₜₒv2₊P
    @test SII.getname(EIndex(:a=>:b, :P)) == :aₜₒb₊P
    @test SII.getname(EIndex(:a=>2, :P)) == :aₜₒv2₊P
    @test SII.getname(EIndex(1=>:b, :P)) == :v1ₜₒb₊P
end

@testset "test comparioson of NWState/NWPara" begin
    # Create test network
    g = path_graph(2)
    vf = Lib.kuramoto_second()
    ef = Lib.kuramoto_edge()
    nw = Network(g, vf, ef)

    # Test NWParameter comparison
    p1 = NWParameter(nw; default=false)
    pflat(p1) .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    p2 = NWParameter(nw; default=false)
    pflat(p2) .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    p3 = NWParameter(nw; default=false)
    pflat(p3) .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.1]  # slightly different

    # Test exact equality
    @test p1 == p2
    @test p1 != p3

    # Test approximate equality
    @test isapprox(p1, p2)
    @test isapprox(p1, p3; atol=0.2)
    @test !isapprox(p1, p3; atol=0.05)

    # Test NWState comparison
    s1 = NWState(nw; default=false)
    uflat(s1) .= [0.1, 0.2, 0.3, 0.4]  # state values
    pflat(s1) .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]  # parameter values

    s2 = NWState(nw; default=false)
    uflat(s2) .= [0.1, 0.2, 0.3, 0.4]
    pflat(s2) .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    s3 = NWState(nw; default=false)
    uflat(s3) .= [0.1, 0.2, 0.3, 0.41]  # slightly different state
    pflat(s3) .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    # Test exact equality
    @test s1 == s2
    @test s1 != s3

    # Test approximate equality
    @test isapprox(s1, s2)
    @test isapprox(s1, s3; atol=0.02)
    @test !isapprox(s1, s3; atol=0.005)
end

@testset "errormsg quality for indexing" begin
    g = path_graph(2)
    vf = Lib.kuramoto_second()
    ef = Lib.kuramoto_edge()
    nw = Network(g, vf, ef)

    @test nw[VPIndex(1)] === nw[VIndex(1)]
    @test get_default(nw, VPIndex(1, :M)) == get_default(nw, VIndex(1, :M))
end
