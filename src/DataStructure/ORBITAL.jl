
#ORBITAL could be an original data struct constructed from wanner90_centres.xyz.
#If need to define a more elegant struct about wannierbasis, then design a new one.

export ORBITAL, numorb, spinORBITAL
"""
The data in wannier90_centres.xyz.
"""
struct ORBITAL{OT <: Real, AT <: Real}
	location::Vector{CartesianCoordinates{OT}}
	name::Vector{String}
	index::Vector{Int}
	atom_location::Vector{CartesianCoordinates{AT}}
	atom_name::Vector{String}
	belonging::Vector{Int} # The orbitals is belongs to which atom.
end
function ORBITAL(
	location::AbstractVector;
	name = String[],
	index = Int[],
	atom_location = CartesianCoordinates{Float64}[],
	atom_name = String[],
	belonging = Int[],
)

	P = reduce(promote_type, eltype.(location))
	location = map(CartesianCoordinates{P} ∘ collect, location)

	P = reduce(promote_type, eltype.(atom_location))
	atom_location = map(CartesianCoordinates{P} ∘ collect, atom_location)

	num = length(location)

	name = deepcopy(name)
	n = length(name)
	if n ≠ 0 && n ≠ num
		error("Wrong orbital's name.")
	elseif eltype(name) <: AbstractString
		name = string.(name)
	end

	index = deepcopy(index)
	if isempty(index)
		index = collect(1:num)
	elseif length(index) ≠ num
		error("Wrong index of POSCAR.")
	end

	belonging = deepcopy(belonging)
	if isempty(belonging) && !isempty(atom_location)
		belonging = Vector{Int}(undef, num)
		for i in 1:num
			(_, I) = findmin(x -> sum(abs2, x - location[i]), atom_location)
			belonging[i] = I
		end
	end

	return ORBITAL(
		location,
		name,
		index,
		atom_location,
		deepcopy(atom_name),
		belonging,
	)
end

function Base.show(io::IO, orbital::ORBITAL)
	print(io, "ORBITAL with $(numorb(orbital)) orbitals and $(length(orbital.atom_location)) atoms.")
end
numorb(orbital::ORBITAL) = length(orbital.location)


function spinORBITAL(orbital::ORBITAL)::ORBITAL

	belonging = orbital.belonging
	if !isempty(belonging)
		spinbelonging = [belonging; belonging]
	end

	return ORBITAL(
		[orbital.location; orbital.location],
		[orbital.name; orbital.name],
		[orbital.index; orbital.index],
		orbital.atom_location,
		orbital.atom_name,
		spinbelonging,
	)
end
