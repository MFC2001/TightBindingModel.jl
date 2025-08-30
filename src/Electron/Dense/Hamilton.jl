export Hamilton
"""
	Hamilton(k, TB::AbstractTightBindModel)
"""
function Hamilton(k, TB::AbstractTightBindModel)
	n = numorb(TB)
	H = Matrix{ComplexF64}(undef, n, n)

	for j in 1:n, i in 1:j #Access by column.
		H[i, j] = iszero(TB.Nhop[i, j]) ? 0 : sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), TB.hop[i, j])
	end

	return Hermitian(H, :U)
end
"""
	Hamilton(k, TB::AbstractTightBindModel, orblocat::AbstractVector{<:ReducedCoordinates{<:Real}})
"""
function Hamilton(k, TB::AbstractTightBindModel, orblocat::AbstractVector{<:ReducedCoordinates{<:Real}})
	n = numorb(TB)
	H = Matrix{ComplexF64}(undef, n, n)

	for j in 1:n, i in 1:j #Access by column.
		H[i, j] = iszero(TB.Nhop[i, j]) ? 0 : cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), TB.hop[i, j])
	end

	return Hermitian(H, :U)
end
"""
	Hamilton(k, hr::HR)
"""
function Hamilton(k, hr::HR)
	norb = numorb(hr)

	H = Matrix{ComplexF64}(undef, norb, norb)

	for j in 1:norb, i in 1:j
		H[i, j] = iszero(hr.Nhop[i, j]) ? 0 : sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), hr.hop[i, j])
	end

	return Hermitian(H, :U)
end
"""
	Hamilton(k, hr::HR, orblocat::AbstractVector{<:ReducedCoordinates{<:Real}})
"""
function Hamilton(k, hr::HR, orblocat::AbstractVector{<:ReducedCoordinates{<:Real}})
	norb = numorb(hr)

	H = Matrix{ComplexF64}(undef, norb, norb)

	for j in 1:norb, i in 1:j
		H[i, j] = iszero(hr.Nhop[i, j]) ? 0 : cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), hr.hop[i, j])
	end

	return Hermitian(H, :U)
end
function Hamilton(hr::HR)
	norb = numorb(hr)

	H = Matrix{ComplexF64}(undef, norb, norb)

	for j in 1:norb, i in 1:j
		H[i, j] = iszero(hr.Nhop[i, j]) ? 0 : sum(value, hr.hop[i, j])
	end

	return Hermitian(H, :U)
end
