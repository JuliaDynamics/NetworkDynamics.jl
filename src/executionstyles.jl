"""
    abstract type ExecutionStyle{buffered::Bool} end

Abstract type for execution style. The coreloop dispatches based on the
Execution style stored in the network object.

- `buffered=true` means that the edge input es explicitly gathered, i.e. the
  vertex outputs in the output buffer will be copied into a dedicated input buffer
  for the edges.
- `buffered=false` means, that the edge inputs are not explicitly gathered, but
  the corloop will perform a redirected lookup into the output buffer.
"""
abstract type ExecutionStyle{buffered} end

"""
    struct SequentialExecution{buffered::Bool}

Sequential execution, no parallelism. For `buffered` see [`ExecutionStyle`](@ref).
"""
struct SequentialExecution{buffered} <: ExecutionStyle{buffered} end

"""
    struct KAExecution{buffered}

Parallel execution using [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl).
Works with GPU and CPU arrays. For `buffered` see [`ExecutionStyle`](@ref).
"""
struct KAExecution{buffered} <: ExecutionStyle{buffered} end

"""
    struct PolyesterExecution{buffered}

Parallel execution using [`Polyester.jl`](https://github.com/JuliaSIMD/Polyester.jl).
For `buffered` see [`ExecutionStyle`](@ref).
"""
struct PolyesterExecution{buffered} <: ExecutionStyle{buffered} end

"""
    struct ThreadedExecution{buffered}

Parallel execution using Julia threads.
For `buffered` see [`ExecutionStyle`](@ref).
"""
@kwdef struct ThreadedExecution{buffered} <: ExecutionStyle{buffered}
    chunk_cache::Dict{UInt, Vector} = Dict{UInt, Vector}()
end

usebuffer(::ExecutionStyle{buffered}) where {buffered} = buffered
usebuffer(::Type{<:ExecutionStyle{buffered}}) where {buffered} = buffered

# check cuda compatibliity
iscudacompatible(x) = iscudacompatible(typeof(x))
iscudacompatible(::Type{<:Any}) = false
iscudacompatible(::Type{<:KAExecution}) = true
