function BSE_wannier_pp(qgrid, bse, nnkpts, BSEband;
	outfolder = "./",
	bandindex = nothing,
	guess = nothing,
	guessbasis = :iR,
)

	if isnothing(bandindex)
		bandindex = collect(1:length(BSEband[1].values))
	elseif bandindex isa AbstractVector{<:Integer}
		bandindex = Int.(collect(bandindex))
	else
		error("bandindex should be a Vector of Int.")
	end

	M = _mmn_BSE(qgrid, bse, BSEband, bandindex, nnkpts)
	Writemmn(nnkpts, M, joinpath(outfolder, "exciton.mmn"))
	Writeeig(BSEband, joinpath(outfolder, "exciton.eig"); bandindex)

	if !isnothing(guess)
		if guessbasis == :iR
			#require guess(ei,Rₑ,hi,Rₕ)
			W = _guess_iR_to_vck(qgrid, bse, guess)
		elseif guessbasis == :ik
			#TODO
		elseif guessbasis == :vck
			#TODO
		else
			error("Wrong wannierbasis.")
		end
		A = _amn_BSE(BSEband, bandindex, W)
		Writeamn(A, joinpath(outfolder, "exciton.amn"))
	end

	return nothing
end
