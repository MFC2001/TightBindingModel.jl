export WriteBAND
"""
	WriteBAND(data::AbstractVector{LinearAlgebra.Eigen}, file::AbstractString; mode="w", comment="From LatticeModel.WriteBAND.")
	WriteBAND(band::AbstractMatrix{<:Real}, file::AbstractString; mode="w", comment="From LatticeModel.WriteBAND.")
	WriteBAND(kline::KData, data::AbstractVector{LinearAlgebra.Eigen}, file::AbstractString; mode="w", comment="From LatticeModel.WriteBAND.")
	WriteBAND(kline::KData, band::AbstractMatrix{<:Real}, file::AbstractString; mode="w", comment="From LatticeModel.WriteBAND.")

	or for cluster can use WriteBAND(data::Eigen, file::AbstractString; mode="w", comment="From LatticeModel.WriteBAND.")

	size(band) = (Nband, Nk).
"""
function WriteBAND(band::AbstractVector{<:Eigen}, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteBAND.")
	Nk = length(band)
	Nband = length(band[1].values)

	bandmatrix = Matrix{Float64}(undef, Nband, Nk)

	for i in eachindex(band)
		bandmatrix[:, i] = band[i].values
	end

	WriteBAND(bandmatrix, file; mode, comment)
end
function WriteBAND(band::Eigen, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteBAND.")
	WriteBAND(band.values, file; mode, comment)
end
function WriteBAND(band::AbstractVecOrMat{<:Real}, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteBAND.")

	Nk = size(band, 2)
	kline = (; line = 1:Nk, index = [1], name = ["Γ"])

	WriteBAND(kline, band, file; mode, comment)
end
function WriteBAND(kline, band::AbstractVector{<:Eigen}, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteBAND.")
	Nband = length(band[1].values)

	bandmatrix = Matrix{Float64}(undef, Nband, length(kline.line))
	for i in eachindex(band)
		bandmatrix[:, i] = band[i].values
	end

	WriteBAND(kline, bandmatrix, file; mode, comment)
end
function WriteBAND(kline, band::AbstractVecOrMat{<:Real}, file::AbstractString; mode = "w", comment = "From LatticeModel.WriteBAND.")

	Nband = size(band, 1)
	Nk = length(kline.line)

	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = "#" * comment * " " * now

	path = dirname(file)
	mkpath(path)
	file = open(file, mode)

	println(file, comment * " #K-Path(1/A) Energy-Level(eV)")
	println(file, "# NKPTS & NBANDS: ", Nk, " ", Nband)

	write(file, "# kindex: ")
	for index in kline.index
		@printf(file, "%8u", index)
	end
	write(file, "\n# kpoint: ")
	for index in kline.index
		@printf(file, "%8.4f", kline.line[index])
	end
	write(file, "\n# kname:  ")
	for name in kline.name
		@printf(file, "%8s", name)
	end
	println(file, "")

	if Nk == 1
		write(file, "# Only Γ.\n")
		output = Vector{String}(undef, Nband)
		for i in 1:Nband
			output[i] = @sprintf("%8u%14.6f", i, band[i])
		end
		writedlm(file, output)
	else
		line = Vector{String}(undef, Nk)
		for i in 1:Nk
			line[i] = @sprintf("%10.5f", kline.line[i])
		end
		bandi = Vector{String}(undef, Nk)
		for i in 1:Nband
			for j in 1:Nk
				bandi[j] = @sprintf("%14.6f", band[i, j])
			end
			write(file, "# Band-Index  $(i)\n")
			writedlm(file, line .* bandi)
			write(file, "\n")
		end
	end

	close(file)

	return nothing
end
