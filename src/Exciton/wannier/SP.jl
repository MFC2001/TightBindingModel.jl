function BSE_SP_wannier(qgrid::RedKgrid, bse::AbstractBSE;
	nnkpfile::Union{AbstractString, Nothing} = nothing,
	outfolder = "./",
	bandindex = nothing,
	guess = nothing,
	guessbasis = :iR,
)

	if isnothing(nnkpfile)
		nnkpts = kgrid_shell_automatic(qgrid, reciprocal(bse.TB.lattice))
	else
		#TODO
		error("To be continued!")
		nnkpts = Readnnkp(nnkpfile)
	end

	bsefolder = joinpath(outfolder, "BSE")
	mkpath(bsefolder)
	serialize(joinpath(bsefolder, "qgrid.sl"), qgrid)
	serialize(joinpath(bsefolder, "bse.sl"), bse)
	serialize(joinpath(bsefolder, "nnkpts.sl"), nnkpts)

	BSEband = BAND(qgrid, bse; vector = true)
	serialize(joinpath(bsefolder, "bseband.sl"), BSEband)

	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband;
		outfolder = joinpath(outfolder, "triplet"), bandindex, guess, guessbasis)

	return nothing
end
