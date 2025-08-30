
export Cell, atomnames, numatom


abstract type AbstractCell end

@struct_hash_equal_isequal struct Cell{P <: Union{CartesianCoordinates{<:Real}, ReducedCoordinates{<:Real}}} <: AbstractCell
	lattice::Lattice{Float64}
	location::Vector{P}
	name::Vector{String}
	index::Vector{Int}
	period::Vec3{Bool}
end

"""
	Cell(lattice, location; location_type, name = String[], index = Int[], period = Bool[1, 1, 1])

Create a new cell.

Argument `lattice` is a [`Lattice`](@ref) type.
Fractional atomic locations `location` are given by a vector of ``N`` vectors, where ``N`` is the number of atoms.

Argument `location_type` is \"Reduced\" or \"Cartesian\".
Argument `name` is a list of ``N`` values, where the same kind of name need to be a String.
Argument `index` is a list of ``N`` values, where the same kind of index need to be a Int.
Argument `period` = [1(0),1(0),1(0)]

Make sure the basis at unperiodic direction is perpendicular to other basises.
"""
function Cell(lattice, location::AbstractVector; location_type, name = String[], index = Int[], period = [1, 1, 1])

	num = length(location)

	n = length(name)
	if n ≠ 0 && n ≠ num
		throw(DimensionMismatch("The lengths of atomic locations and atomic names are different!"))
	end
	name = string.(name)

	index = deepcopy(index)
	if isempty(index)
		index = collect(1:num)
	elseif length(index) ≠ num
		throw(DimensionMismatch("The lengths of atomic locations and atomic indexs are different!"))
	end

	lattice = deepcopy(lattice)
	if !(lattice isa Lattice)
		lattice = Lattice(lattice)
	end


	P = reduce(promote_type, eltype.(location))
	if location_type[1] ∈ ['R', 'r']
		P = ReducedCoordinates{P}
	elseif location_type[1] ∈ ['C', 'c']
		P = CartesianCoordinates{P}
	else
		error("`location_type` only can be set as \"Reduced\" or \"Cartesian\".")
	end
	location = map(P ∘ collect, location)

	period = Vec3{Bool}(collect(period))

	return Cell{P}(lattice, location, name, index, period)
end

atomnames(cell::Cell) = unique(cell.name)
numatom(cell::Cell) = length(cell.location)
Lattice(cell::Cell) = deepcopy(cell.lattice)

function Base.sort(cell::Cell{P}) where {P}
	T = sortslices([cell.name cell.location collect(1:length(cell.name))]; dims = 1)

	name = similar(cell.name)
	name .= T[:, 1]

	location = similar(cell.location)
	location .= T[:, 2]

	index = cell.index[T[:, 3]]

	return Cell{P}(deepcopy(cell.lattice), location, name, index, deepcopy(cell.period))
end

function Base.convert(cell::Cell{ReducedCoordinates{P}}) where {P <: Real}

	location = map(x -> cell.lattice * x, cell.location)

	return Cell{CartesianCoordinates{P}}(deepcopy(cell.lattice), location, deepcopy(cell.name), deepcopy(cell.index), deepcopy(cell.period))
end
function Base.convert(cell::Cell{CartesianCoordinates{P}}) where {P <: Real}

	location = map(x -> cell.lattice \ x, cell.location)

	return Cell{ReducedCoordinates{P}}(deepcopy(cell.lattice), location, deepcopy(cell.name), deepcopy(cell.index), deepcopy(cell.period))
end
