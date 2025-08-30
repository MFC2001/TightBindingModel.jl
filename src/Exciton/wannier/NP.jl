function BSE_NP_wannier(qgrid::MonkhorstPack, bse::AbstractBSE; kwards...)
	return BSE_NP_wannier(RedKgrid(qgrid), bse; kwards...)
end
function BSE_NP_wannier(qgrid::RedKgrid, bse::AbstractBSE;
	nnkpfile::Union{AbstractString, Nothing} = nothing,
	outfolder = "./",
	bandindex_t = nothing,
	bandindex_s = nothing,
	guess_t = nothing,
	guess_s = nothing,
	guessbasis_t = :iR,
	guessbasis_s = :iR,
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

	(BSEband_t, BSEband_s) = BAND(qgrid, bse; vector = true)
	serialize(joinpath(bsefolder, "triplet.sl"), BSEband_t)
	serialize(joinpath(bsefolder, "singlet.sl"), BSEband_s)

	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband_t;
		outfolder = joinpath(outfolder, "triplet"), bandindex = bandindex_t, guess = guess_t, guessbasis = guessbasis_t)
	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband_s;
		outfolder = joinpath(outfolder, "singlet"), bandindex = bandindex_s, guess = guess_s, guessbasis = guessbasis_s)

	return nothing
end
