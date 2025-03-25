module NetworkDynamics
using Graphs: Graphs, AbstractGraph, SimpleEdge, edges, vertices, ne, nv,
              SimpleGraph, SimpleDiGraph, add_edge!, has_edge
using TimerOutputs: TimerOutputs, @timeit_debug, @timeit

using ArgCheck: @argcheck
using PreallocationTools: PreallocationTools, DiffCache, get_tmp
using SciMLBase: SciMLBase
using Base.Threads: @threads
using NNlib: NNlib
using KernelAbstractions: KernelAbstractions, @kernel, @index, @Const, get_backend
using Atomix: Atomix
using Polyester: Polyester
using Mixers: Mixers
using LinearAlgebra: LinearAlgebra, UniformScaling
using SparseArrays: sparse
using StyledStrings: StyledStrings, @styled_str
using RecursiveArrayTools: RecursiveArrayTools, DiffEqArray
using FastClosures: @closure
using ForwardDiff: ForwardDiff
using Printf: @sprintf
using Random: Random
using Static: Static, StaticInt
using SciMLBase: VectorContinuousCallback, CallbackSet, DiscreteCallback
using DiffEqCallbacks: DiffEqCallbacks

@static if VERSION ≥ v"1.11.0-0"
    using Base: AnnotatedIOBuffer, AnnotatedString
else
    using StyledStrings: AnnotatedIOBuffer, AnnotatedString
end

using Base: @propagate_inbounds
using InteractiveUtils: subtypes

import SymbolicIndexingInterface as SII
using StaticArrays: StaticArrays, SVector

include("utils.jl")

export VertexModel, EdgeModel
export StateMask, Symmetric, AntiSymmetric, Directed, Fiducial
export FeedForwardType, PureFeedForward, FeedForward, NoFeedForward, PureStateMap
export dim, sym, pdim, psym, obssym, hasinsym, insym, hasindim, indim,
       outdim, outsym, fftype, metadata, symmetadata
include("component_functions.jl")

export ExecutionStyle, SequentialExecution, KAExecution, ThreadedExecution, PolyesterExecution
include("executionstyles.jl")

export Network
include("network_structure.jl")

export Aggregator, KAAggregator, SequentialAggregator, PolyesterAggregator, ThreadedAggregator, SparseAggregator
export ff_to_constraint
include("aggregators.jl")
include("gbufs.jl")
include("construction.jl")
include("coreloop.jl")
include("external_inputs.jl")

# XXX: have both, s[:] and uflat(s) ?
export VIndex, EIndex, VPIndex, EPIndex, NWState, NWParameter, uflat, pflat
export vidxs, eidxs, vpidxs, epidxs
export save_parameters!
export @obsex
include("symbolicindexing.jl")

export has_metadata, get_metadata, set_metadata!
export has_default, get_default, set_default!, set_defaults!
export has_guess, get_guess, set_guess!
export has_init, get_init, set_init!
export has_bounds, get_bounds, set_bounds!
export has_graphelement, get_graphelement, set_graphelement!
export get_initial_state, dump_initial_state, dump_state
export has_callback, get_callbacks, set_callback!, add_callback!
export has_position, get_position, set_position!
export has_marker, get_marker, set_marker!
include("metadata.jl")

export ComponentCondition, ComponentAffect
export ContinousComponentCallback, VectorContinousComponentCallback
export DiscreteComponentCallback, PresetTimeComponentCallback
export SymbolicView
include("callbacks.jl")

using NonlinearSolve: AbstractNonlinearSolveAlgorithm, NonlinearFunction
using NonlinearSolve: NonlinearLeastSquaresProblem, NonlinearProblem
using SteadyStateDiffEq: SteadyStateProblem, SteadyStateDiffEqAlgorithm, SSRootfind
export find_fixpoint, initialize_component!, init_residual, set_interface_defaults!
include("initialization.jl")

include("show.jl")

const CHECK_COMPONENT = Ref(true)
export chk_component
include("doctor.jl")

#=
using StyledStrings
s1 = styled"{bright_red:brred} {bright_green:brgreen} {bright_yellow:bryellow} {bright_blue:brblue} {bright_magenta:brmagenta} {bright_cyan:brcyan} {bright_black:brblack} {bright_white:brwhite}";
s2 = styled"{red:  red} {green:  green} {yellow:  yellow} {blue:  blue} {magenta:  magenta} {cyan:  cyan} {black:  black} {white:  white}";
println(s1, "\n", s2)
=#
const ND_FACES = [
    :NetworkDynamics_inactive => StyledStrings.Face(foreground=:bright_black),
    :NetworkDynamics_defaultval => StyledStrings.Face(weight=:light),
    :NetworkDynamics_guessval => StyledStrings.Face(foreground=:bright_black, weight=:light),
    :NetworkDynamics_fordstsrc => StyledStrings.Face(foreground=:bright_blue),
    :NetworkDynamics_fordst => StyledStrings.Face(foreground=:bright_yellow),
    :NetworkDynamics_forsrc => StyledStrings.Face(foreground=:bright_magenta),
    :NetworkDynamics_forlayer => StyledStrings.Face(foreground=:bright_blue),
    :NetworkDynamics_name => StyledStrings.Face(weight=:bold),
    :NetworkDynamics_fftype => StyledStrings.Face(foreground=:bright_blue),
]

__init__() = foreach(StyledStrings.addface!, ND_FACES)

function reloadfaces!()
    StyledStrings.resetfaces!()
    for (k,v) in NetworkDynamics.ND_FACES
        delete!(StyledStrings.FACES.default, k)
    end
    foreach(StyledStrings.addface!, NetworkDynamics.ND_FACES)
end
# NetworkDynamics.reloadfaces!()

using PrecompileTools: @compile_workload
@compile_workload begin
    include("precompile_workload.jl")
end

end
