function checkKeywords(finduc, outparameter)
	if finduc ∉ ["auto", "translate", "custom"]
		error("Wrong keyword finduc.")
	elseif outparameter[1] ∉ ["value", "index"]
		error("Wrong keyword outhr.")
	elseif outparameter[2] ∉ ["value", "index"]
		error("Wrong keyword outorb.")
	elseif outparameter[3][1] ∉ ['Y', 'y', 'N', 'n']
		error("Wrong keyword outucorientation.")
	end
	return nothing
end

function checkCell(cell1, cell2)
	if Set(cell1.name) ≠ Set(cell2.name)
		error("Please check if unitcell and supercell are compatible.")
	end

	location_type = typeof(cell1).parameters[1]
	if location_type <: ReducedCoordinates{<:Real}
		cell1 = convert(cell1)
	elseif location_type <: CartesianCoordinates{<:Real}
		cell1 = deepcopy(cell1)
	end

	location_type = typeof(cell2).parameters[1]
	if location_type <: ReducedCoordinates{<:Real}
		cell2 = convert(cell2)
	elseif location_type <: CartesianCoordinates{<:Real}
		cell2 = deepcopy(cell2)
	end

	return cell1, cell2
end
