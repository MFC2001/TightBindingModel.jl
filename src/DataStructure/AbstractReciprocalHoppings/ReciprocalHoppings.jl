
abstract type ReciprocalHoppings{T <: Number} <: AbstractReciprocalHoppings end

@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates)
	A = Matrix{ComplexF64}(undef, rh.norb, rh.norb)
	return rh(A, k)
end
@inline function (rh::ReciprocalHoppings)(k::ReducedCoordinates, orblocat)
	A = Matrix{ComplexF64}(undef, rh.norb, rh.norb)
	return rh(A, k, orblocat)
end
@inline function (rh::ReciprocalHoppings)(A::AbstractMatrix, parameters...)
	size(A) == size(rh.hops) || error("Buffer size mismatch.")
	map(CartesianIndices(rh.hops)) do I
		A[I] = rh(I, parameters...)
	end
	return A
end
@inline function (rh::ReciprocalHoppings)(I::CartesianIndex{2}, k::ReducedCoordinates)
	iszero(rh.Nhop[I]) ? 0 :
	sum(hop -> hop.t * hopphase(hop, k), rh.hops[I])
end
@inline function (rh::ReciprocalHoppings)(I::CartesianIndex{2}, k::ReducedCoordinates, orblocat)
	(i, j) = Tuple(I)
	iszero(rh.Nhop[i, j]) ? 0 :
	cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), rh.hops[i, j])
end
@inline function (rh::ReciprocalHoppings)(i::Integer, j::Integer, k::ReducedCoordinates)
	iszero(rh.Nhop[i, j]) ? 0 :
	sum(hop -> hop.t * hopphase(hop, k), rh.hops[i, j])
end
@inline function (rh::ReciprocalHoppings)(i::Integer, j::Integer, k::ReducedCoordinates, orblocat)
	iszero(rh.Nhop[i, j]) ? 0 :
	cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * hopphase(hop, k), rh.hops[i, j])
end

Base.zero(::U, norb::Integer) where {T <: Number, U <: ReciprocalHoppings{T}} = zero(U, norb)
function Base.zero(::Type{U}, norb::Integer) where {T <: Number, U <: ReciprocalHoppings{T}}
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	Nhop = zeros(Int, norb, norb)
	return U(norb, hops, Nhop)
end

function HR(rh::ReciprocalHoppings{T}) where {T}
	N = sum(rh.Nhop)
	path = Matrix{Int}(undef, 5, N)
	value = Vector{T}(undef, N)

	n = 0
	for I in CartesianIndices(rh.hops)
		for hop in rh.hops[I]
			n += 1
			path[:, n] = [hop.R; hop.i; hop.j]
			value[n] = hop.t
		end
	end

	return HR{T}(collect(1:rh.norb), transpose(path), value)
end

function _buildhops_from_hr(norb, path::AbstractMatrix, value::AbstractVector{T}) where {T}
	hops = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	for i in axes(path, 1)
		push!(hops[path[i, 4], path[i, 5]], Hopping(path[i, :], value[i]))
	end
	return hops
end


