
function symmetry_operations(TB::TightBindModel; symtol = 1e-5, check_symmetry = true)
	return symmetry_operations(TB.basis, TB.atom_name, TB.atom_location; symtol, check_symmetry)
end
function symmetry_operations(cell::Cell; symtol = 1e-5, check_symmetry = true)
	return symmetry_operations(parent(cell.lattice), cell.name, cell.location; symtol, check_symmetry)
end
