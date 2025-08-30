export vckMap

struct vckMap
	nv::Int
	nc::Int
	nk::Int
	v2idx::Dict{Int, Int}
	c2idx::Dict{Int, Int}
	idx2v::Vector{Int}
	idx2c::Vector{Int}
	vck2idx::Array{Int, 3}
	idx2vck::Vector{Tuple{Int, Int, Int}}
end

function vckMap(v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, nk::Integer)

	idx2v = map(Int, v)
	idx2c = map(Int, c)

	nv = length(idx2v)
	nc = length(idx2c)
	nk = Int(nk)

	v2idx = Dict{Int, Int}()
	for i in eachindex(idx2v)
		v2idx[idx2v[i]] = i
	end
	c2idx = Dict{Int, Int}()
	for i in eachindex(idx2c)
		c2idx[idx2c[i]] = i
	end

	idx2vck = Vector{Tuple{Int, Int, Int}}(undef, nv * nc * nk)
	vck2idx = Array{Int, 3}(undef, nv, nc, nk)
	n = 0
	for k in 1:nk, c in 1:nc, v in 1:nv
		n += 1
		idx2vck[n] = (idx2v[v], idx2c[c], k)
		vck2idx[v, c, k] = n
	end

	return vckMap(nv, nc, nk, v2idx, c2idx, idx2v, idx2c, vck2idx, idx2vck)
end
vckMap(v::AbstractVector{<:Integer}, c::Integer, nk::Integer) = vckMap(v, [c], nk)
vckMap(v::Integer, c::AbstractVector{<:Integer}, nk::Integer) = vckMap([v], c, nk)
vckMap(v::Integer, c::Integer, nk::Integer) = vckMap([v], [c], nk)

Base.length(vckmap::vckMap) = length(vckmap.idx2vck)
Base.size(vckmap::vckMap) = (vckmap.nv, vckmap.nc, vckmap.nk)
function Base.getindex(vckmap::vckMap, idx::Integer)
	return vckmap.idx2vck[idx]
end
function Base.getindex(vckmap::vckMap, vck::Tuple{<:Integer, <:Integer, <:Integer})
	return vckmap.vck2idx[vckmap.v2idx[vck[1]], vckmap.c2idx[vck[2]], vck[3]]
end


# struct vckVector{T <: Number} <: AbstractVector{T}
# 	data::Array{T, 3}
# 	vckmap::vckMap
# end
# struct vckMatrix{T <: Number} <: AbstractMatrix{T}
# 	data::Array{T, 6}
# 	vckmap::vckMap
# end


# function vckVector(
# 	v_vals::AbstractVector{<:Integer},
# 	c_vals::AbstractVector{<:Integer},
# 	nk::Integer,
# 	data = Array{T, 3}(undef, length(v_vals), length(c_vals), nk),
# )

# end

# Base.size(H::BSEHamiltonian) = (H.nv, H.nc, H.nk, H.nv, H.nc, H.nk)
# Base.getindex(H::HermitianStorage, i::Int, j::Int) = begin
# 	if H.part == :U
# 		i ≤ j ? H.storage[triu_index(i, j, H.n)] : conj(H.storage[triu_index(j, i, H.n)])
# 	else
# 		i ≥ j ? H.storage[tril_index(i, j, H.n)] : conj(H.storage[tril_index(j, i, H.n)])
# 	end
# end
