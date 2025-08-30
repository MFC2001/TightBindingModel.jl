
function WSfactor_grid(lattice::Lattice, kgrid::MonkhorstPack, realgrid)

	#划分程度，取奇数
	N₁ = realgrid[1]
	N₂ = realgrid[2]
	N₃ = realgrid[3]
	N = N₁ * N₂ * N₃

	n₁ = Int((N₁ - 1) // 2)
	n₂ = Int((N₂ - 1) // 2)
	n₃ = Int((N₃ - 1) // 2)


	N_half = Int((N - 1) // 2)
	points_frac = Vector{SVector{3, Rational{Int}}}(undef, N_half)
	n = 0
	for i in 1:n₁, j in -n₂:n₂, k in -n₃:n₃
		n += 1
		points_frac[n] = SVector{3}(i // N₁, j // N₂, k // N₃)
	end
	for j in 1:n₂, k in -n₃:n₃
		n += 1
		points_frac[n] = SVector{3}(0 // N₁, j // N₂, k // N₃)
	end
	for k in 1:n₃
		n += 1
		points_frac[n] = SVector{3}(0 // N₁, 0 // N₂, k // N₃)
	end

	# r = map(r -> 1 / norm(lattice * r), points_frac)
	inv_r = map(r -> 1 / norm(lattice * r), points_frac)

	Nk = length(kgrid)
	Ω = det(lattice.data)
	dV = Nk * Ω / N
	rsz = (3 * dV / 4π)^(1 / 3)

	I_factor = dV / 4π

	rlattice = reciprocal(lattice)

	I = function (k)
		k2 = sum(x -> x^2, rlattice * k)
		k1 = sqrt(k2)
		return sum(i -> inv_r[i] * cos(2π * (k ⋅ points_frac[i])), eachindex(inv_r)) * k2 * I_factor * 2 +
			   1 - cos(k1 * rsz)
	end

end
