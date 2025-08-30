
export MonkhorstPack, RedKgrid, IrredKgrid, Kline

abstract type AbstractBrillouinZone end

"""
Perform BZ sampling employing a Monkhorst-Pack grid.
"""
struct MonkhorstPack <: AbstractBrillouinZone
	kgrid_size::Vec3{Int}
	kshift::Vec3{Rational{Int}}
	function MonkhorstPack(size, kshift)
		map(kshift) do ks
			ks in (0, 1 // 2) || error("Only kshifts of 0 or 1//2 implemented.")
		end
		new(size, kshift)
	end
end
MonkhorstPack(kgrid_size::AbstractVector; kshift = [0, 0, 0]) = MonkhorstPack(kgrid_size, kshift)
MonkhorstPack(kgrid_size::Tuple; kshift = [0, 0, 0]) = MonkhorstPack(kgrid_size, kshift)
MonkhorstPack(k1::Integer, k2::Integer, k3::Integer) = MonkhorstPack([k1, k2, k3])
function Base.show(io::IO, kgrid::MonkhorstPack)
	print(io, "MonkhorstPack(", kgrid.kgrid_size)
	if !iszero(kgrid.kshift)
		print(io, ", ", Float64.(kgrid.kshift))
	end
	print(io, ")")
end
Base.length(kgrid::MonkhorstPack) = prod(kgrid.kgrid_size)


"""Bring ``k``-point coordinates into the range [-0.5, 0.5)"""
function normalize_kdirect(x::Real)
	x = x - round(x, RoundNearestTiesUp)
	@assert -0.5 ≤ x < 0.5
	return x
end
normalize_kdirect(k::AbstractVector{<:Real}) = normalize_kdirect.(k)


""" """
struct RedKgrid <: AbstractBrillouinZone
	kdirect::Vector{ReducedCoordinates{Rational{Int}}}
	kgrid_size::Vec3{Int}
	kshift::Vec3{Rational{Int}}
end
function RedKgrid(kgrid::MonkhorstPack)
	start = -floor.(Int, (kgrid.kgrid_size .- 1) .// 2)
	stop = ceil.(Int, (kgrid.kgrid_size .- 1) .// 2)

	kgrid_index = [-kgrid.kshift .+ Vec3([i, j, k]) for i ∈ start[1]:stop[1], j ∈ start[2]:stop[2], k ∈ start[3]:stop[3]]
	kgrid_index = reshape(kgrid_index, :)
	kdirect = [x .// kgrid.kgrid_size for x in kgrid_index]

	return RedKgrid(normalize_kdirect.(kdirect), kgrid.kgrid_size, kgrid.kshift)
end
function Base.show(io::IO, redkgrid::RedKgrid)
	print(io, "Reducibal kgrids with $(length(redkgrid.kdirect)) reducible k-points.")
end

Base.getindex(redkgrid::RedKgrid, index...) = getindex(redkgrid.kdirect, index...)
# Base.setindex!(redkgrid::RedKgrid, v, i::Int) = (redkgrid.kdirect[i] = v)
# Base.size(redkgrid::RedKgrid) = size(redkgrid.kdirect)
Base.length(redkgrid::RedKgrid) = length(redkgrid.kdirect)

# Base.axes(redkgrid::RedKgrid) = axes(redkgrid.kdirect)

Base.iterate(redkgrid::RedKgrid, state = 1) = state > length(redkgrid) ? nothing : (redkgrid[state], state + 1)
Base.eltype(::RedKgrid) = ReducedCoordinates{Rational{Int}}

Base.firstindex(::RedKgrid) = 1
Base.lastindex(redkgrid::RedKgrid) = length(redkgrid)
Base.eachindex(redkgrid::RedKgrid) = eachindex(redkgrid.kdirect)

# Base.view(redkgrid::RedKgrid, inds...) = view(redkgrid.kdirect, inds...)

Base.keys(redkgrid::RedKgrid) = keys(redkgrid.kdirect)




struct IrredKgrid{UT <: Number} <: AbstractBrillouinZone
	kdirect::Vector{ReducedCoordinates{Rational{Int}}}
	kweight::Vector{Rational{Int}}
	redkdirect::Vector{ReducedCoordinates{Rational{Int}}}
	irmap::Vector{Int}
	symop::Vector{SymOp}
	lattsymop::Vector{LattSymOp{UT}}
	kgrid_size::Vec3{Int}
	kshift::Vec3{Rational{Int}}
end
function Base.show(io::IO, irredkgrid::IrredKgrid)
	print(io, "Irreducibal kgrids with $(length(irredkgrid.kdirect)) irreducible k-points",
		" and $(length(irredkgrid.redkdirect)) reducible k-points.")
end
Base.getindex(irredkgrid::IrredKgrid, index...) = getindex(irredkgrid.kdirect, index...)
Base.length(irredkgrid::IrredKgrid) = length(irredkgrid.kdirect)
Base.iterate(irredkgrid::IrredKgrid, state = 1) = state > length(irredkgrid) ? nothing : (irredkgrid[state], state + 1)
Base.eltype(::IrredKgrid) = ReducedCoordinates{Rational{Int}}
Base.firstindex(::IrredKgrid) = 1
Base.lastindex(irredkgrid::IrredKgrid) = length(irredkgrid)
Base.eachindex(irredkgrid::IrredKgrid) = eachindex(irredkgrid.kdirect)
Base.keys(irredkgrid::IrredKgrid) = keys(irredkgrid.kdirect)


"""
Explicitly define the k-points along which to perform BZ sampling.
(Useful for bandstructure calculations)
"""
struct ExplicitKpoints{T} <: AbstractBrillouinZone
	kcoords::Vector{Vec3{T}}
	kweights::Vector{T}
end
function ExplicitKpoints(kcoords::AbstractVector{<:AbstractVector{T}}, kweights::AbstractVector{T}) where {T}
	@assert length(kcoords) == length(kweights)
	ExplicitKpoints{T}(kcoords, kweights)
end
function ExplicitKpoints(kcoords::AbstractVector{<:AbstractVector{T}}) where {T}
	ExplicitKpoints(kcoords, ones(T, length(kcoords)) ./ length(kcoords))
end
function Base.show(io::IO, kgrid::ExplicitKpoints)
	print(io, "ExplicitKpoints with $(length(kgrid.kcoords)) k-points")
end
Base.length(kgrid::ExplicitKpoints) = length(kgrid.kcoords)




"""
Kline includes kpoints and names.
"""
struct Kline{L <: Real} <: AbstractBrillouinZone
	rlattice::ReciprocalLattice{L}
	kdirect::Vector{ReducedCoordinates{Float64}}
	line::Vector{Float64}
	name::Vector{String}
	index::Vector{Int}
end
function Kline(; basis = I, kdirect::AbstractVector, line = Float64[], name = String[], index = Int[])

	try
		if length(kdirect) ≠ length(line) && !any(isempty, [kdirect, line])
			throw(DimensionMismatch("The lengths of kcoords and klines are different!"))
		end
	catch err
		println(err.msg, " Try to set kline = 1:Nk.")
		line = Float64.(collect(1:length(kdirect)))
	end

	P = reduce(promote_type, eltype.(kdirect))
	P = ReducedCoordinates{P}
	kdirect = map(P ∘ collect, kdirect)

	if length(name) ≠ length(index)
		throw(DimensionMismatch("The lengths of high symmetry points and its index are different!"))
	end

	name = string.(name)

	return Kline(ReciprocalLattice(basis), kdirect, line, name, index)
end
Base.getindex(kline::Kline, index...) = getindex(kline.kdirect, index...)
Base.length(kline::Kline) = length(kline.kdirect)
Base.iterate(kline::Kline, state = 1) = state > length(kline) ? nothing : (kline[state], state + 1)
Base.eltype(::Kline) = ReducedCoordinates{Float64}
Base.firstindex(::Kline) = 1
Base.lastindex(kline::Kline) = length(kline)
Base.eachindex(kline::Kline) = eachindex(kline.kdirect)
Base.keys(kline::Kline) = keys(kline.kdirect)
