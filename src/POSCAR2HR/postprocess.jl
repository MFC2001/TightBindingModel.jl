function postprocessHR(outhr, hr_path, datahr, orbindex)
	if outhr == "value"

		hr_value = datahr.value[hr_path[:, 6]]
		hr_path = hr_path[:, 1:5]

		hr = HR(hr_path, hr_value; orbindex, hrsort = "Y")

	elseif outhr == "index"
		hr = hr_path
	end
	return hr
end
function postprocessORBITAL(outorb, orbital_index, dataorbital, datalattice, rotatelattice, aimcell)
	#orbital_index:[atom_index atom_orb_index data_atom_index data_orb_index]
	if outorb == "value"
		orbital_dlocation_frac = map(i -> datalattice \ (dataorbital.location[i] - dataorbital.atom_location[dataorbital.belonging[i]]), eachindex(dataorbital.location))

		orbital_name = aimcell.name[orbital_index[:, 1]]
		orbital_location = map((ai, uoi) -> aimcell.location[ai] + rotatelattice * orbital_dlocation_frac[uoi], orbital_index[:, 1], orbital_index[:, 4])

		orbital = ORBITAL(orbital_location; name = orbital_name, index = orbital_index[:, 4],
			atom_location = aimcell.location, atom_name = aimcell.name, belonging = orbital_index[:, 1])

	elseif outorb == "index"
		orbital = orbital_index
	end
	return orbital
end
