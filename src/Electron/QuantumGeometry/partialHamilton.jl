
function partialHamilton(k, lattice::Lattice, hr::HR, orblocat::AbstractVector)
	norb = numorb(hr)

	H = Array{ComplexF64}(undef, norb, norb, 3)

	for j in 1:norb, i in 1:j
		dorb = orblocat[j] - orblocat[i]
		H[i, j, :] = iszero(hr.Nhop[i, j]) ? zeros(ComplexF64, 3) : lattice * (im * cis(2π * (k ⋅ dorb)) * sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)) * (hop.R + dorb), hr.hop[i, j]))
		H[j, i, :] = conj.(H[i, j, :])
	end

	return H
end
