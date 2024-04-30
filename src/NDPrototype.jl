module NDPrototype
using Graphs
using OrderedCollections
using Unrolled: unrolled_foreach, @unroll
using TimerOutputs: @timeit_debug, reset_timer!, print_timer

using ArgCheck: @argcheck
using PreallocationTools: LazyBufferCache
using EnumX: @enumx
using SciMLBase: SciMLBase
using Base.Threads: @threads
using NNlib: NNlib
using KernelAbstractions: KernelAbstractions, @kernel, @index, @Const, get_backend
using Atomix: Atomix

include("utils.jl")
include("edge_coloring.jl")
include("component_functions.jl")
export Network, SequentialExecution, ThreadedExecution
include("network_structure.jl")

export NaiveAggregator, NNlibScatter, KAAggregator, SequentialAggregator
include("aggregators.jl")
include("construction.jl")
include("coreloop.jl")

end
