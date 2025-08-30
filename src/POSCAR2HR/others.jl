
function classifyCell(cell::Cell, unitcell::Cell, deps::Real)
	#=
	Find atom index for every kind of atoms in cell. 
	unitcell's every atom is its own kind, because function `namelist!` have changed unitcell.name.
	=#
	deps = deps^2

	supercell = SuperCell(unitcell)

	name = similar(cell.name)
	for i in eachindex(cell.location)

		aimlocation = cell.location[i]

		D = map(x -> norm2(x - aimlocation), supercell.location)
		(minD, aim) = findmin(D)

		if minD > deps

			presentN = length(supercell.location)
			ExpandSuperCell!(supercell, aimlocation)
			D = map(x -> norm2(x - aimlocation), supercell.location[presentN+1:end])

			(minD, aim) = findmin(D)
			aim = aim + presentN

			if minD > deps
				@show √minD
				error("The cell and unitcell may be incompatible.")
			end
		end
		name[i] = supercell.name[aim]
	end


	atom_index = similar(unitcell.name, Vector{Int})
	for i in eachindex(unitcell.name)
		atom_index[i] = findall(name .== unitcell.name[i])
	end

	return atom_index
end


function namelist!(cell::Cell)
	#Will reset the name and index of cell.
	elem_name = unique(cell.name)
	elem_num = zeros(Int, length(elem_name))

	for i in eachindex(cell.name)
		I = findfirst(x -> x == cell.name[i], elem_name)
		elem_num[I] += 1
		cell.name[i] = cell.name[i] * string(elem_num[I])
		cell.index[i] = elem_num[I]
	end

	return nothing
end

function norm2(v::AbstractVector{<:Number})::Real
	return v ⋅ v
end

function calc_vecangle(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real})::Real
	return atan(norm(v1 × v2), v1 ⋅ v2)
end
