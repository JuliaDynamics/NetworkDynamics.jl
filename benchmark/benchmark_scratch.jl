using Dates
using NetworkDynamics
using Chairmarks
using Graphs
using Serialization
using StableRNGs
using Test
using Random
using GLMakie
using JET
using SciMLBase
(isinteractive() ? includet : include)("benchmark_utils.jl")
(isinteractive() ? includet : include)("benchmark_models.jl")

#####################################
# serialization and deserialization #
#####################################
serialize(Dates.format(now(), "yyyy-mm-dd")*"_results.data", bd)

# bd = deserialize("./2024-07-08_results.data")
bd = deserialize("./2024-06-25_results.dat")

fig = plot_over_N(bd)
DataInspector(fig)
comp = bd


#############################
# single call for profiling #
#############################
N = 10
g = watts_strogatz(N, N ÷ 2, 0.0; rng=StableRNG(1))
edge = diffusion_edge()
# edge = diffusion_dedge()
vertex = diffusion_vertex()
# execution = PolyesterExecution{false}()
execution = KAExecution{false}()
# execution = SequentialExecution{true}()
# aggregator = PolyesterAggregator(+)
aggregator = SequentialAggregator(+)
_nd = Network(g, vertex, edge; execution, aggregator)
_x0 = rand(dim(_nd))
_dx = similar(_x0)

@report_opt ignored_modules=(AnyFrameModule(NetworkDynamics.PreallocationTools),) _nd(_dx, _x0, nothing, NaN)
@report_call ignored_modules=(AnyFrameModule(NetworkDynamics.PreallocationTools),)  _nd(_dx, _x0, nothing, NaN)

@descend _nd(_dx, _x0, nothing, NaN)
@b $_nd($_dx,$_x0, nothing, NaN)

####

N = 10
vert, edg, g = heterogeneous(N)
execution = KAExecution{false}()
# execution = ThreadedExecution{false}()
aggregator = PolyesterAggregator(+)
_nd = Network(g, vert, edg; execution, aggregator)
_x0 = rand(dim(_nd));
_dx = similar(_x0);
_p = rand(pdim(_nd));

@report_opt ignored_modules=(AnyFrameModule(NetworkDynamics.PreallocationTools),) _nd(_dx, _x0, _p, NaN)
@report_call ignored_modules=(AnyFrameModule(NetworkDynamics.PreallocationTools),)  _nd(_dx, _x0, _p, NaN)

@time _nd(_dx, _x0, _p, NaN)
@b $_nd($_dx,$_x0, $_p, NaN)

### cuda
using NetworkDynamics
using CUDA
using NetworkDynamics.Adapt
N = 10
vert, edg, g = heterogeneous(N)
execution = KAExecution{true}()
aggregator = KAAggregator(+)
_nd = Network(g, vert, edg; execution, aggregator)
_x0 = rand(dim(_nd));
_dx = similar(_x0);
_p = rand(pdim(_nd));
_nd(_dx, _x0, _p, NaN)

_x0_d = CuArray(_x0)
_dx_d  = CuArray(_dx)
_p_d = CuArray(_p)

_nd_d = adapt(_x0_d, _nd)

_nd(_dx_d, _x0_d, _p_d, NaN)


Vector(_dx_d) ≈ _dx

isbits(Type{Float64})
isbits(Float64)
isbits(1.0)

@b $_nd($_dx,$_x0, $_p, NaN)



####
function evalnd(nd::T, dx, x0) where {T}
    h = hash(dx)
    for i in 1:10
        nd(dx, x0, nothing, NaN)
        h = hash(dx, h)
    end
    h
end
@b evalnd($_nd,$_dx,$_x0)

Profile.clear()
@pprof evalnd(_nd, _dx, _x0)

@profview evalnd(_nd, _dx, _x0)

@pprof @b $_nd($_dx, $_x0, nothing, NaN)

@descend _nd(_dx, _x0, nothing, NaN)
