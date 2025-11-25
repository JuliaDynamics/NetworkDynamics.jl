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
using SparseArrays: SparseArrays, sparse, SparseMatrixCSC
using StyledStrings: StyledStrings, @styled_str
using RecursiveArrayTools: RecursiveArrayTools, DiffEqArray
using FastClosures: @closure
using ForwardDiff: ForwardDiff
using Printf: @sprintf
using Random: Random
using Static: Static, StaticInt
using SciMLBase: VectorContinuousCallback, CallbackSet, DiscreteCallback
using DiffEqCallbacks: DiffEqCallbacks
using MacroTools: postwalk, @capture

@static if VERSION ≥ v"1.11.0-0"
    using Base: AnnotatedIOBuffer, AnnotatedString
else
    using StyledStrings: AnnotatedIOBuffer, AnnotatedString
end

using Base: @propagate_inbounds
using InteractiveUtils: subtypes

import SymbolicIndexingInterface as SII
using SymbolicIndexingInterface: variable_symbols, parameter_symbols
using StaticArrays: StaticArrays, SVector

export implicit_output, ComponentPostprocessing
include("utils.jl")

export VertexModel, EdgeModel
export StateMask, Symmetric, AntiSymmetric, Directed, Fiducial
export FeedForwardType, PureFeedForward, FeedForward, NoFeedForward, PureStateMap
export dim, sym, pdim, psym, obssym, hasinsym, insym, hasindim, indim,
       outdim, outsym, fftype, metadata, symmetadata
include("component_functions.jl")

export ExecutionStyle, SequentialExecution, KAExecution, ThreadedExecution, PolyesterExecution
include("executionstyles.jl")

export Network, get_graph, set_jac_prototype!
include("network_structure.jl")

export Aggregator, KAAggregator, SequentialAggregator, PolyesterAggregator, ThreadedAggregator, SparseAggregator
export ff_to_constraint
include("aggregators.jl")
include("gbufs.jl")
include("construction.jl")
include("coreloop.jl")
include("external_inputs.jl")

# XXX: have both, s[:] and uflat(s) ?
export VIndex, EIndex, VPIndex, EPIndex
export ParamIdx, StateIdx
export save_parameters!, extract_nw
export variable_symbols, parameter_symbols
include("symbolicindexing_base.jl")

export NWState, NWParameter, uflat, pflat
export vidxs, eidxs, vpidxs, epidxs, generate_indices, FilteringProxy
export @obsex
include("symbolicindexing.jl")

export ComponentCondition, ComponentAffect
export ContinuousComponentCallback, VectorContinuousComponentCallback
export DiscreteComponentCallback, PresetTimeComponentCallback
export SymbolicView
include("callbacks.jl")

export @initconstraint, InitConstraint
export @initformula, InitFormula
export @guessformula, GuessFormula
include("init_constraints.jl")

using OrderedCollections: OrderedDict
using NonlinearSolve: NonlinearFunction, NonlinearLeastSquaresProblem
using SteadyStateDiffEq: SteadyStateProblem, SSRootfind
export find_fixpoint, set_interface_defaults!
export initialize_component, initialize_component!, init_residual
export initialize_componentwise, initialize_componentwise!, interface_values
export NetworkInitError, ComponentInitError
include("initialization.jl")

export has_metadata, get_metadata, set_metadata!, delete_metadata!, strip_metadata!
export has_default, get_default, set_default!, delete_default!, set_defaults!, strip_defaults!
export has_guess, get_guess, set_guess!, delete_guess!, strip_guesses!
export has_init, get_init, set_init!, delete_init!, strip_inits!
export has_bounds, get_bounds, set_bounds!, delete_bounds!, strip_bounds!
export has_graphelement, get_graphelement, set_graphelement!
export get_initial_state, dump_initial_state, dump_state
export has_callback, get_callbacks, set_callback!, add_callback!, delete_callbacks!
export has_initconstraint, get_initconstraints, set_initconstraint!, add_initconstraint!, delete_initconstraints!
export has_initformula, get_initformulas, set_initformula!, add_initformula!, delete_initformulas!
export has_guessformula, get_guessformulas, set_guessformula!, add_guessformula!, delete_guessformulas!
export has_position, get_position, set_position!
export has_marker, get_marker, set_marker!
export get_defaults_dict, get_guesses_dict, get_bounds_dict, get_inits_dict
export free_p, free_u
include("metadata.jl")

export isfixpoint, is_linear_stable, jacobian_eigenvals
include("linear_stability.jl")

include("show.jl")

const CHECK_COMPONENT = Ref(true)
export chk_component
include("doctor.jl")

# additional utils, which depend on specific types beeing defined so they should be
# loaded after all other files
include("post_utils.jl")

export describe_vertices, describe_edges, get_jac_prototype
function describe_vertices end
function describe_edges end
function get_jac_prototype end
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

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f ∈ (describe_vertices, describe_edges)
                ext = Base.get_extension(NetworkDynamics, :NetworkDynamicsDataFramesExt)
                if isnothing(ext)
                    printstyled(io, "\nLoad `DataFrames` in order to use `describe_vertices` and `describe_edges`."; bold=true, color=:red)
                end
            elseif exc.f ∈ (get_jac_prototype,)
                ext = Base.get_extension(NetworkDynamics, :NetworkDynamicsSparsityExt)
                if isnothing(ext)
                    printstyled(io, "\nLoad `SparseConnectivityTracer` in order to use `get_jac_prototype`."; bold=true, color=:red)
                end
            end
        end
    end

    foreach(StyledStrings.addface!, ND_FACES)
end

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

####
#### Deprecated functionality (remove in breaking version)
####
export ContinousComponentCallback, VectorContinousComponentCallback
include("deprecated.jl")

end
