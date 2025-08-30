export WriteHR
"""
	WriteHR(hr::HR, file::AbstractString; mode="w", comment="From MyWrite.Hr.")

Write HR to file.
"""
function WriteHR(hr::HR, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteHR.")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)
	write(file, comment * "\n")

	@printf(file, "%12u\n", numorb(hr))
	allpath = union(x[1:3] for x in eachrow(hr.path))
	nkpoint = length(allpath)
	@printf(file, "%12u\n", nkpoint)

	for i in 1:nkpoint
		@printf(file, "%5u", 1)
		if mod(i, 15) == 0
			@printf(file, "\n")
		end
	end
	if mod(nkpoint, 15) ≠ 0
		@printf(file, "\n")
	end

	for i in eachindex(hr.value)
		@printf(file, "%5u %5u %5u %7u %7u %12.6f %12.6f\n", hr.path[i, :]..., reim(hr.value[i])...)
	end

	close(file)
	return nothing
end

WriteHR(TB::AbstractTightBindModel, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteHR.") =
	WriteHR(HR(TB), file; mode, comment)
