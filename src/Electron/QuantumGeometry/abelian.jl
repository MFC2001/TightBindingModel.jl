export QuantumGeometry
function QuantumGeometry(kgrid::AbstractVector, hr::HR, lattice::Lattice, orblocat::AbstractVector, n::Integer = 1)

	norb = numorb(hr)
	Nk = length(kgrid)

	QG = Array{ComplexF64}(undef, Nk, 3, 3)

	sumindex = setdiff(1:norb, n)

	for k in 1:Nk
		kpoint = kgrid[k]
		E = eigen!(Hamilton(kpoint, hr, orblocat))
		pH = partialHamilton(kpoint, lattice, hr, orblocat)

		E₀ = E.values[n]
		φ₀ = E.vectors[:, n]
		for μ in 1:3, ν in 1:3
			QG[k, μ, ν] = sum(sumindex) do index
				return (φ₀ ⋅ (pH[:, :, μ] * E.vectors[:, index])) * (E.vectors[:, index] ⋅ (pH[:, :, ν] * φ₀)) / (E₀ - E.values[index])^2
			end
		end
	end

	return QG
end
