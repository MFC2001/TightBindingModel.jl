export WritePOSCAR
"""
	WritePOSCAR(poscar::POSCAR, file::AbstractString; mode="w", format="Cartesian", comment="From LatticeModel.WritePOSCAR.")

Write POSCAR to file.
"""
function WritePOSCAR(cell::Cell, file::AbstractString; mode = "w", format = "Cartesian", comment = "From LatticeModel.WritePOSCAR.")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)
	write(file, comment * "\n")
	write(file, "1.0\n")

	for i in 1:3
		@printf(file, "%26.16f %26.16f %26.16f\n", cell.lattice[:, i]...)
	end

	elem_name = atomnames(cell)
	for name in elem_name
		@printf(file, " %5s", name)
	end
	write(file, "\n")

	for name in elem_name
		@printf(file, " %5u", count(cell.name .== name))
	end
	write(file, "\n")

	if format[1] ∈ ['C', 'c']
		write(file, "Cartesian\n")
		for r in cell.location
			@printf(file, "%23.16f %23.16f %23.16f\n", cell.lattice * r...)
		end
	elseif format[1] ∈ ['D', 'd']
		write(file, "Direct\n")
		for r in cell.location
			@printf(file, "%23.16f %23.16f %23.16f\n", r...)
		end
	end

	close(file)

	return nothing
end

WritePOSCAR(TB::AbstractTightBindModel, file::AbstractString; mode = "w", format = "Cartesian", comment = "From LatticeModel.WritePOSCAR.") =
	WritePOSCAR(Cell(TB), file; mode, format, comment)
