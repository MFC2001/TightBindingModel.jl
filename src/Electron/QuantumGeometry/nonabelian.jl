export QuantumGeometry
function QuantumGeometry(kgrid::AbstractVector, hr::HR, lattice::Lattice, orblocat::AbstractVector, n::AbstractVector{<:Integer} = [1])

	norb = numorb(hr)
	nn = length(n)
	Nk = length(kgrid)

	QG = Array{ComplexF64}(undef, Nk, 3, 3, nn, nn)

	sumindex = setdiff(1:norb, n)

	Threads.@threads for k in 1:Nk
		kpoint = kgrid[k]
		E = eigen!(Hamilton(kpoint, hr, orblocat))
		pH = partialHamilton(kpoint, lattice, hr, orblocat)


		for μ in 1:3, ν in 1:3, i in 1:nn, j in 1:nn
			ni = n[i]
			nj = n[j]
			QG[k, μ, ν, i, j] = sum(sumindex) do index
				return (E.vectors[:, ni] ⋅ (pH[:, :, μ] * E.vectors[:, index])) * (E.vectors[:, index] ⋅ (pH[:, :, ν] * E.vectors[:, nj])) /
					   ((E.values[ni] - E.values[index]) * (E.values[nj] - E.values[index]))
			end
		end
	end

	return QG
end
