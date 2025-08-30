export DensityMatrix_k, DensityMatrix_R
function DensityMatrix_R(redkgrid::RedKgrid, Distribution::Function)

	Nk = length(redkgrid)
	R = gridindex(redkgrid.kgrid_size)
	# minusmap = kgridmap(redkgrid, -)

	DM = function (band_red)

		norb = length(band_red[1].values)

		ρk = Array{ComplexF64}(undef, norb, norb, Nk)
		Threads.@threads for k in 1:Nk
			fvector = Distribution(band_red[k].values)

			wavek = band_red[k].vectors
			fwave = wavek .* transpose(fvector)

			# ρk[:, :, k] = conj.(wavek) * transpose(fwave)

			for i in axes(wavek, 1), j in 1:i
				ρk[i, j, k] = wavek[i, :] ⋅ fwave[j, :]
				ρk[j, i, k] = conj(ρk[i, j, k])
			end
		end

		ρR = Array{ComplexF64}(undef, norb, norb, Nk)
		for (Ri, Rv) in enumerate(R)
			for i in 1:norb, j in 1:norb
				# ρR[i, j, Ri] = sum((ki, kv) -> ρk[i, j, ki] * cis(2π * (kv ⋅ Rv)), enumerate(redkgrid.kdirect)) / Nk
				ρR[i, j, Ri] = sum(ki -> ρk[i, j, ki] * cis(2π * (redkgrid.kdirect[ki] ⋅ Rv)), 1:Nk) / Nk
			end
		end

		# ρ = function (i, R₁, j, R₂)
		# 	return ρR[i, j, minusmap[R₂, R₁]]
		# end

		return ρR
	end

	return DM
end
function DensityMatrix_R(irredkgrid::IrredKgrid, Distribution::Function)

	Nk = length(irredkgrid.redkdirect)
	R = gridindex(redkgrid.kgrid_size)

	DM = function (band_irred)
		norb = length(band_irred[1].values)

		ρk = Array{ComplexF64}(undef, norb, norb, Nk)
		Threads.@threads for k in 1:Nk
			fvector = Distribution(band_irred[irredkgrid.irmap[k]].values)

			wavek = irredkgrid.lattsymop[k].U * band_irred[irredkgrid.irmap[k]].vectors
			fwave = wavek .* transpose(fvector)

			for i in axes(wavek, 1), j in 1:i
				ρk[i, j, k] = wavek[i, :] ⋅ fwave[j, :]
				ρk[j, i, k] = conj(ρk[i, j, k])
			end
		end

		ρR = Array{ComplexF64}(undef, norb, norb, Nk)
		Threads.@threads for (Ri, Rv) in enumerate(R)
			for i in 1:norb, j in 1:norb
				ρR[i, j, Ri] = sum((ki, kv) -> ρk[i, j, ki] * cis(2π * (kv ⋅ Rv)), enumerate(irredkgrid.redkdirect)) / Nk
			end
		end

		return ρR
	end

	return DM
end



function DensityMatrix_k(irredkgrid::IrredKgrid, Distribution::Function)

	DM = function (band_irred)
		norb = length(band_irred[1].values)

		ρ = Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}(undef, length(irredkgrid.redkdirect))
		Threads.@threads for k in eachindex(ρ)
			fvector = Distribution(band_irred[irredkgrid.irmap[k]].values)

			wavek = irredkgrid.lattsymop[k].U * band_irred[irredkgrid.irmap[k]].vectors
			fwave = wavek .* transpose(fvector)

			ρk = Matrix{ComplexF64}(undef, norb, norb)
			for j in 1:norb, i in 1:j
				ρk[i, j] = wavek[i, :] ⋅ fwave[j, :]
			end

			ρ[k] = Hermitian(ρk, :U)
		end

		return ρ
	end

	return DM
end
function DensityMatrix_k(::Type{RedKgrid}, Distribution::Function)

	DM = function (band_red)
		norb = length(band_red[1].values)

		ρ = Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}(undef, length(band_red))
		Threads.@threads for k in eachindex(ρ)
			fvector = Distribution(band_red[k].values)

			wavek = band_red[k].vectors
			fwave = wavek .* transpose(fvector)

			ρk = Matrix{ComplexF64}(undef, norb, norb)
			for j in 1:norb, i in 1:j
				ρk[i, j] = wavek[i, :] ⋅ fwave[j, :]
			end

			ρ[k] = Hermitian(ρk, :U)
		end

		return ρ
	end

	return DM
end

function DensityMatrix_k!(irredkgrid::IrredKgrid, Distribution::Function)

	DM! = function (ρ, band)
		norb = length(band[1].values)

		Threads.@threads for k in eachindex(ρ)
			fvector = Distribution(band[irredkgrid.irmap[k]].values)

			wavek = irredkgrid.lattsymop[k].U * band[irredkgrid.irmap[k]].vectors
			fwave = wavek .* transpose(fvector)

			ρk = Matrix{ComplexF64}(undef, norb, norb)
			for j in 1:norb, i in 1:j
				ρk[i, j] = wavek[i, :] ⋅ fwave[j, :]
			end

			ρ[k] = Hermitian(ρk, :U)
		end

		return ρ
	end

	return DM!
end
