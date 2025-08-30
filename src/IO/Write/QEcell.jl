export WriteQEcell
"""
	WritePOSCAR(poscar::POSCAR, file::AbstractString; mode="w", format="Cartesian", comment="From LatticeModel.WriteQEcell.")

Write POSCAR to file.
"""
function WriteQEcell(cell::Cell, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteQEcell.")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)
	write(file, comment * "\n")
	write(file, "\n")


	write(file, "CELL_PARAMETERS (angstrom)\n")
	for i in 1:3
		@printf(file, "%26.16f %26.16f %26.16f\n", cell.lattice[:, i]...)
	end


	write(file, "ATOMIC_POSITIONS (crystal)\n")
	for i in eachindex(cell.name)
		@printf(file, "%-6s %20.10f %20.10f %20.10f\n", cell.name[i], cell.location[i]...)
	end


	close(file)

	return nothing
end
