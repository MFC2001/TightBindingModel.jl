
function Writeeig(band::AbstractVector{<:Eigen}, file::AbstractString;
	bandindex::Union{Nothing, AbstractVector{<:Integer}} = nothing, mode = "w", comment = "From LatticeModel.Writeeig.")

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	if isnothing(bandindex)
		bandindex = eachindex(band[1].values)
	end


	path = dirname(file)
	mkpath(path)
	file = open(file, mode)


	for (i, E) in enumerate(band)
		for (ii, ei) in enumerate(bandindex)
			@printf(file, "%12u %12u %22.12f \n", ii, i, E.values[ei])
		end
	end

	close(file)

	return nothing
end
