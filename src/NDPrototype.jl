module NDPrototype
using Graphs
using OrderedCollections
using Unrolled

using ArgCheck: @argcheck
using PreallocationTools: LazyBufferCache
using EnumX: @enumx
using SciMLBase: SciMLBase
using Base.Threads: @threads
using NNlib: NNlib

include("utils.jl")
include("edge_coloring.jl")
include("component_functions.jl")
include("network_structure.jl")

export NaiveAggregator, NNlibScatter
include("aggregators.jl")
include("construction.jl")
include("coreloop.jl")

end
