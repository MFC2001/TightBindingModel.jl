export WriteORBITAL
"""
	WriteORBITAL(orbital::ORBITAL, file::AbstractString; mode="w", comment="From LatticeModel.WriteORBITAL.")

	WriteORBITAL(orbital::ORBITAL, poscar::POSCAR, file::AbstractString; mode="w", comment="From LatticeModel.WriteORBITAL.")

Write ORBITAL to file.
"""
function WriteORBITAL(orbital::ORBITAL, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteORBITAL.")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)

	N = string(numorb(orbital) + length(orbital.atom_location))
	write(file, N, "\n")
	write(file, comment * "\n")

	if isempty(orbital.index)
		for i in eachindex(orbital.location)
			@printf(file, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "X", orbital.location[i]...)
		end
	else
		for i in eachindex(orbital.location)
			@printf(file, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "X$(orbital.index[i])", orbital.location[i]...)
		end
	end

	try
		!isempty(orbital.atom_location) || error("The atom_location od orbital is empty.")
		for i in eachindex(orbital.atom_location)
			@printf(file, "%-5s\t%15.8f\t%15.8f\t%15.8f\n", "$(orbital.atom_name[i])", orbital.atom_location[i]...)
		end
	catch e
		println(file, e.msg)
	end

	close(file)

	return nothing
end

WriteORBITAL(TB::AbstractTightBindModel, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteORBITAL.") =
	WriteORBITAL(ORBITAL(TB), file; mode, comment)
