module DFTK_irrk

using StaticArrays
using LinearAlgebra

const Mat3{T} = SMatrix{3, 3, T, 9} where {T}
const Vec3{T} = SVector{3, T} where {T}



include("./bzmesh.jl")




# function time_reversal(spg_mesh)
# 	kmesh = spg_mesh.mesh
# 	kcoords = map(x -> [x[1] // kmesh[1], x[2] // kmesh[2], x[3] // kmesh[3]], spg_mesh.grid_address)
# 	kcoords_normalized = normalize_kpoint_coordinate.(kcoords)

# 	rspg_mesh = deepcopy(spg_mesh)

# 	for i in eachindex(kcoords_normalized)
# 		@show aim = findindex(normalize_kpoint_coordinate(-kcoords_normalized[i]), kcoords_normalized)
# 		if rspg_mesh.ir_mapping_table[i] ≠ rspg_mesh.ir_mapping_table[aim]
# 			I = rspg_mesh.ir_mapping_table .== rspg_mesh.ir_mapping_table[aim]
# 			rspg_mesh.ir_mapping_table[I] .= rspg_mesh.ir_mapping_table[i]
# 		end
# 	end

# 	return rspg_mesh
# end

# findindex(x::AbstractArray{<:Rational}, X)  = findall(isequal(x), X)
# findindex(x::AbstractArray{T}, X) where {T} = findall(y -> isapprox(x, y; atol = sqrt(eps(T))), X)

end
