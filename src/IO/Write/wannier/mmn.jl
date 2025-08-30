
function Writemmn(nnkpts::AbstractMatrix{<:Integer}, M::AbstractArray{<:Number}, file::AbstractString; mode = "w", comment = "From LatticeModel.Writemmn.")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)
	write(file, comment * "\n")

	nband = size(M, 1)
	Nk = maximum(nnkpts[1, :])
	nb = Int(size(nnkpts, 2) // Nk)

	@printf(file, "%12u %12u %12u \n", nband, Nk, nb)

	for i in axes(M, 3)
		@printf(file, "%5u %5u %5u %5u %5u \n", nnkpts[:, i]...)
		for n in Base.OneTo(nband)
			for m in Base.OneTo(nband)
				@printf(file, "%18.13f %18.13f \n", reim(M[m, n, i])...)
			end
		end
	end

	close(file)

	return nothing
end
