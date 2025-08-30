function vck2iR(bse::BSE, q::AbstractVector)

	Nk = length(bse.kgrid)
	bandk = BAND(bse.kgrid, bse.TB; vector = true)
	bandkq = BAND(map(k -> k + q, bse.kgrid.kdirect), bse.TB; vector = true)

	btM = function (vi, ci, ki, ei, Rₑ, hi, Rₕ)
		return bandkq[ki].vectors[ei, bse.c[ci]] * conj(bandk[ki].vectors[hi, bse.v[vi]]) * cis((kgrid.kdirect[ki] + q) * Rₑ - kgrid.kdirect[ki] * Rₕ)
	end
end

function iR2vck()

end
