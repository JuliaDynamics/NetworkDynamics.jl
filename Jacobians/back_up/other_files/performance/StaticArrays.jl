using LinearAlgebra
using StaticArrays

SA[1, 2, 3] isa SVector{3,Int}

v1 = SVector(1, 2, 3)

#struct JacGraphData1
#    v_jac_array::Array{Array{Float64, 2}, 1}
#    e_jac_array::Array{Array{Array{Float64, 2}, 1}, 1}
#    e_jac_product::Array{Array{Float64, 1}, 1}

v_dims = 3
v_jac_array = [Array{Float64,2}(undef, dim, dim) for dim in v_dims]
a = @SArray randn(2, 2, 2, 2, 2, 2)
zeros(SVector{3})
@SArray randn(2)

SArray{Tuple{3,3},Int64,2,9}
SArray{Tuple{3,3},Int64,9}
SArray{Tuple{3,3}}

SArray()
SArray()
