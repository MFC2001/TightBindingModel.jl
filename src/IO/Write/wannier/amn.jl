
function Writeamn(A, file::AbstractString; mode = "w", comment = "From LatticeModel.Writeamn.")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	(nband, nw, Nk) = size(A)

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)
	write(file, comment * "\n")

	@printf(file, "%12u %12u %12u \n", nband, Nk, nw)

	for k in 1:Nk
		for wi in 1:nw
			for n in 1:nband
				@printf(file, "%5u %5u %5u %18.12f %18.12f \n", n,wi,k,reim(A[n, wi, k])...)
			end
		end
	end

	close(file)

	return nothing
end
