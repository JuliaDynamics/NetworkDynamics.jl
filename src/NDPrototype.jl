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
using Polyester: Polyester
using Mixers: Mixers
using Parameters: @with_kw_noshow
using LinearAlgebra: LinearAlgebra
using DocStringExtensions
using StyledStrings: StyledStrings, @styled_str


using Adapt: Adapt, adapt

include("utils.jl")
include("edge_coloring.jl")
include("component_functions.jl")
export Network, SequentialExecution, KAExecution
include("network_structure.jl")

export NaiveAggregator, NNlibScatter, KAAggregator, SequentialAggregator,
       PolyesterAggregator
include("aggregators.jl")
include("construction.jl")
include("coreloop.jl")

include("adapt.jl")
include("show.jl")

#=
styled"{bright_red:red} {bright_green:green} {bright_yellow:yellow} {bright_blue:blue} {bright_magenta:magenta} {bright_cyan:cyan}"
styled"{red:red} {green:green} {yellow:yellow} {blue:blue} {magenta:magenta} {cyan:cyan}"
=#
const ND_FACES = [
    :NetworkDynamics_defaultval => StyledStrings.Face(foreground=:bright_black),
    :NetworkDynamics_fordstsrc => StyledStrings.Face(foreground=:bright_blue),
    :NetworkDynamics_fordst => StyledStrings.Face(foreground=:bright_yellow),
    :NetworkDynamics_forsrc => StyledStrings.Face(foreground=:bright_magenta),
    :NetworkDynamics_forlayer => StyledStrings.Face(foreground=:bright_blue),
]

__init__() = foreach(StyledStrings.addface!, ND_FACES)

function reloadfaces!()
    StyledStrings.resetfaces!()
    for (k,v) in NDPrototype.ND_FACES
        delete!(StyledStrings.FACES.default, k)
    end
    foreach(StyledStrings.addface!, NDPrototype.ND_FACES)
end

end
