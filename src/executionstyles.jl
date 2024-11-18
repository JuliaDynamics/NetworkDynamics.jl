abstract type ExecutionStyle{buffered} end
struct SequentialExecution{buffered} <: ExecutionStyle{buffered} end
struct KAExecution{buffered} <: ExecutionStyle{buffered} end
struct PolyesterExecution{buffered} <: ExecutionStyle{buffered} end
struct ThreadedExecution{buffered} <: ExecutionStyle{buffered} end
usebuffer(::ExecutionStyle{buffered}) where {buffered} = buffered
usebuffer(::Type{<:ExecutionStyle{buffered}}) where {buffered} = buffered

# check cuda compatibliity
iscudacompatible(x) = iscudacompatible(typeof(x))
iscudacompatible(::Type{<:Any}) = false
iscudacompatible(::Type{<:KAExecution}) = true
