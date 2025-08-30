
abstract type AbstractTightBindModel <: AbstractModel end

#TODO 如何将自旋指标引入。

struct TightBindModel{T <: Number, U <: ReciprocalHoppings{T}} <: AbstractTightBindModel
	lattice::Lattice{Float64}
	atom_name::Vector{String}
	atom_location::Vector{ReducedCoordinates{Float64}}
	orb_name::Vector{String}
	orb_location::Vector{ReducedCoordinates{Float64}}
	H::HermitianReciprocalHoppings{T, U}
	period::Vec3{Bool}
end

function TightBindModel(hr::HR, cell::Cell{<:ReducedCoordinates}, orbital::ORBITAL)
	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orbital.name,
		[cell.lattice \ x for x in orbital.location],
		HermitianReciprocalHoppings(hr),
		cell.period,
	)
end
function TightBindModel(hr::HR, cell::Cell{<:CartesianCoordinates}, orbital::ORBITAL)
	return TightBindModel(hr, convert(cell), orbital)
end

function Base.show(io::IO, TB::TightBindModel)
	print(io, "$(count(TB.period)) dimensinal Tight binding model with $(numatom(TB)) atoms and $(numorb(TB)) orbitals.")
end
numatom(TB::TightBindModel) = length(TB.atom_location)
numorb(TB::TightBindModel) = length(TB.orb_location)

# Base.convert(::Type{TightBindModel{T₁}}, TB::TightBindModel{T₂}) where {T₁, T₂} = 
# 	TightBindModel{T₁}(TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location,
# 		convert(HermitianReciprocalHoppings{T₁}, TB.H), TB.period)


Cell(TB::TightBindModel) = Cell{ReducedCoordinates{Float64}}(TB.lattice, TB.atom_location, TB.atom_name, eachindex(TB.atom_location), TB.period)
HR(TB::TightBindModel) = HR(TB.H)
function ORBITAL(TB::AbstractTightBindModel)
	return ORBITAL(
		map(x -> TB.lattice * x, TB.orb_location);
		name = TB.orb_name,
		atom_location = map(x -> TB.lattice * x, TB.atom_location),
		atom_name = TB.atom_name,
	)
end

function spinTightBindModel(TB::TightBindModel{T, U}) where {T <: Number, U <: ReciprocalHoppings{T}}
	spin_orb_name = [TB.orb_name .* 'u'; TB.orb_name .* 'd']
	spin_orb_location = [TB.orb_location; TB.orb_location]

	norb = length(TB.orb_name)
	hops_up = deepcopy(TB.H.rh.hops)
	hops_dn = map(hops_up) do hops
		map(hop -> Hopping{T}(hop.i, hop.j, hop.R, conj(hop.t)), hops)
	end
	hops_0 = [Vector{Hopping{T}}(undef, 0) for _ in CartesianIndices((norb, norb))]
	spinhops = [
		hops_up hops_0;
		hops_0 hops_dn
	]
	spinNhop = length.(spinhops)
	spinrh = BaseReciprocalHoppings{T}(norb * 2, spinNhop, spinhops)
	return TightBindModel(
		TB.lattice,
		TB.atom_name,
		TB.atom_location,
		spin_orb_name,
		spin_orb_location,
		Hermitian(spinrh),
		TB.period,
	)
end
