export coulomb_inner_product
function coulomb_inner_product(irredkgrid::IrredKgrid, Uhr::HR, Uhr_onsite::HR, Jhr::HR)::Function

	U = LatticeCoulomb(Uhr)
	U_onsite = LatticeCoulomb(Uhr_onsite)
	J = LatticeCoulomb(Jhr)

	Nk = length(irredkgrid.redkdirect)
	Uk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Uk_onsite = Vector{Matrix{ComplexF64}}(undef, Nk)
	Jk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Threads.@threads for k in 1:Nk
		Uk[k] = U(irredkgrid.redkdirect[k]) / Nk
		Uk_onsite[k] = U_onsite(irredkgrid.redkdirect[k]) / Nk
		Jk[k] = J(irredkgrid.redkdirect[k]) / Nk
	end

	minusmap = kgridmap(irredkgrid, -)
	addmap = kgridmap(irredkgrid, +)

	Vik(i, k₁, j, k₂, k, k₃, l, k₄; spin = false) =
		spin ?
		(i == l && j == k) * Uk_onsite[minusmap[k₃, k₂]][i, j] +
		(i == k && j == l) * Jk[minusmap[k₄, k₂]][i, j] +
		(i == j && k == l) * Jk[addmap[k₁, k₂]][i, k] :
		(i == l && j == k) * Uk[minusmap[k₃, k₂]][i, j] +
		(i == k && j == l) * Jk[minusmap[k₄, k₂]][i, j] +
		(i == j && k == l) * Jk[addmap[k₁, k₂]][i, k]

	# Vik(i, σ₁, k₁, j, σ₂, k₂, k, σ₃, k₃, l, σ₄, k₄) =
	# (i == l && j == k && σ₁ == σ₄ && σ₂ == σ₃) * Uk[minusmap[k₃, k₂]][i, j] +
	# (i == k && j == l) * Jk[minusmap[k₄, k₂]][i, j] +
	# (i == j && k == l) * Jk[addmap[k₁, k₂]][i, k]

	return Vik
end
function coulomb_inner_product(irredkgrid::IrredKgrid, Uhr::HR, Jhr::HR)::Function

	U = LatticeCoulomb(Uhr)
	J = LatticeCoulomb(Jhr)

	Nk = length(irredkgrid.redkdirect)
	Uk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Jk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Threads.@threads for k in 1:Nk
		Uk[k] = U(irredkgrid.redkdirect[k]) / Nk
		Jk[k] = J(irredkgrid.redkdirect[k]) / Nk
	end

	minusmap = kgridmap(irredkgrid, -)
	addmap = kgridmap(irredkgrid, +)

	Vik(i, k₁, j, k₂, k, k₃, l, k₄) =
		(i == l) * (j == k) * Uk[minusmap[k₃, k₂]][i, j] +
		(i == k) * (j == l) * Jk[minusmap[k₄, k₂]][i, j] +
		(i == j) * (k == l) * Jk[addmap[k₁, k₂]][i, k]

	return Vik
end
function coulomb_inner_product(redkgrid::RedKgrid, Uhr::HR, Jhr::HR)::Function

	U = LatticeCoulomb(Uhr)
	J = LatticeCoulomb(Jhr)

	Nk = length(redkgrid.kdirect)
	Uk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Jk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Threads.@threads for k in Base.OneTo(Nk)
		Uk[k] = U(redkgrid.kdirect[k]) / Nk
		Jk[k] = J(redkgrid.kdirect[k]) / Nk
	end

	minusmap = kgridmap(redkgrid, -)
	addmap = kgridmap(redkgrid, +)


	Vik(i, k₁, j, k₂, k, k₃, l, k₄) =
		(i == l) * (j == k) * Uk[minusmap[k₃, k₂]][i, j] +
		(i == k) * (j == l) * Jk[minusmap[k₄, k₂]][i, j] +
		(i == j) * (k == l) * Jk[addmap[k₁, k₂]][i, k]

	return Vik
end

