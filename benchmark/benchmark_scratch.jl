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

┌ (::Network{…})(du::Vector{…}, u::Vector{…}, p::Nothing, t::Float64) @ NetworkDynamics /home/hw/.julia/packages/TimerOutputs/Lw5SP/src/TimerOutput.jl:253
│┌ process_layer!(::KAExecution{…}, nw::Network{…}, layer::NetworkDynamics.NetworkLayer{…}, dupt::Tuple{…}) @ NetworkDynamics /home/hw/.julia/dev/NetworkDynamics/src/coreloop.jl:155
││┌ unrolled_foreach(f::NetworkDynamics.var"#122#123"{…}, t::Tuple{…}) @ NetworkDynamics /home/hw/.julia/dev/NetworkDynamics/src/utils.jl:17
│││┌ (::NetworkDynamics.var"#122#123"{…})(batch::NetworkDynamics.EdgeBatch{…}) @ NetworkDynamics /home/hw/.julia/dev/NetworkDynamics/src/coreloop.jl:158
││││┌ kwcall(::@NamedTuple{…}, ::KernelAbstractions.Kernel{…}, ::Type{…}, ::@NamedTuple{…}, ::Vector{…}, ::Vector{…}, ::Vector{…}, ::Vector{…}, ::Vector{…}, ::Nothing, ::Float64) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:39
│││││┌ (::KernelAbstractions.Kernel{…})(::Type{…}, ::@NamedTuple{…}, ::Vector{…}, ::Vector{…}, ::Vector{…}, ::Vector{…}, ::Vector{…}, ::Nothing, ::Float64; ndrange::Int64, workgroupsize::Nothing) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:46
││││││┌ __run(obj::KernelAbstractions.Kernel{…}, ndrange::Tuple{…}, iterspace::KernelAbstractions.NDIteration.NDRange{…}, args::Tuple{…}, dynamic::KernelAbstractions.NDIteration.DynamicCheck, static_threads::Bool) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:84
│││││││┌ __thread_run(tid::Int64, len::Any, rem::Any, obj::KernelAbstractions.Kernel{…}, ndrange::Tuple{…}, iterspace::KernelAbstractions.NDIteration.NDRange{…}, args::Tuple{…}, dynamic::KernelAbstractions.NDIteration.DynamicCheck) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:117
││││││││┌ cpu_ekernel!(__ctx__::KernelAbstractions.CompilerMetadata{…} where _A, ::Type, batch::@NamedTuple{…}, du::Vector{…}, u::Vector{…}, s::Vector{…}, srcrange::Vector{…}, dstrange::Vector{…}, p::Nothing, t::Float64) @ NetworkDynamics /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/macros.jl:285
│││││││││┌ __validindex(ctx::KernelAbstractions.CompilerMetadata{…} where _A, idx::CartesianIndex{…}) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:156
││││││││││┌ expand(ndrange::KernelAbstractions.NDIteration.NDRange{…}, groupidx::Integer, idx::CartesianIndex{…}) @ KernelAbstractions.NDIteration /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/nditeration.jl:92
│││││││││││┌ getindex(A::CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, I::Integer) @ Base ./abstractarray.jl:1315
││││││││││││┌ to_indices(A::CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, I::Tuple{Integer}) @ Base ./indices.jl:365
│││││││││││││┌ to_indices(A::CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, inds::Tuple{}, I::Tuple{Integer}) @ Base ./indices.jl:368
││││││││││││││ runtime dispatch detected: Base.to_index(A::CartesianIndices{1, Tuple{Base.OneTo{Int64}}}, %2::Integer)::Any
│││││││││││││└────────────────────
│││││││││┌ __validindex(ctx::KernelAbstractions.CompilerMetadata{…} where _A, idx::CartesianIndex{…}) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:156
││││││││││ runtime dispatch detected: KernelAbstractions.expand(%2::KernelAbstractions.NDIteration.NDRange{…}, %3::Any, idx::CartesianIndex{1})::CartesianIndex{1}
│││││││││└────────────────────
││││││││┌ cpu_ekernel!(__ctx__::KernelAbstractions.CompilerMetadata{…} where _A, ::Type, batch::@NamedTuple{…}, du::Vector{…}, u::Vector{…}, s::Vector{…}, srcrange::Vector{…}, dstrange::Vector{…}, p::Nothing, t::Float64) @ NetworkDynamics /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/macros.jl:286
│││││││││┌ __index_Global_Linear(ctx::KernelAbstractions.CompilerMetadata{…} where _A, idx::CartesianIndex{…}) @ KernelAbstractions /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/cpu.jl:137
││││││││││ runtime dispatch detected: KernelAbstractions.expand(%1::KernelAbstractions.NDIteration.NDRange{…}, %2::Any, idx::CartesianIndex{1})::CartesianIndex{1}
│││││││││└────────────────────
││││││││┌ cpu_ekernel!(__ctx__::KernelAbstractions.CompilerMetadata{…} where _A, ::Type, batch::@NamedTuple{…}, du::Vector{…}, u::Vector{…}, s::Vector{…}, srcrange::Vector{…}, dstrange::Vector{…}, p::Nothing, t::Float64) @ NetworkDynamics /home/hw/.julia/packages/KernelAbstractions/MAxUm/src/macros.jl:287
│││││││││┌ apply_edge_unbuffered!(::Type, batch::@NamedTuple{…}, i::Int64, du::Vector{…}, u::Base.Experimental.Const{…}, s::Vector{…}, srcrange::Base.Experimental.Const{…}, dstrange::Base.Experimental.Const{…}, p::Nothing, t::Float64) @ NetworkDynamics /home/hw/.julia/dev/NetworkDynamics/src/coreloop.jl:174
││││││││││┌ _has_dynamic(T::Type) @ NetworkDynamics /home/hw/.julia/dev/NetworkDynamics/src/coreloop.jl:279
│││││││││││ runtime dispatch detected: NetworkDynamics.statetype(T::Type)::Any
││││││││││└────────────────────
││││││││││┌ _has_dynamic(T::Type) @ NetworkDynamics /home/hw/.julia/dev/NetworkDynamics/src/coreloop.jl:279
│││││││││││ runtime dispatch detected: NetworkDynamics._has_dynamic(%1::Any)::Bool
││││││││││└────────────────────


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
