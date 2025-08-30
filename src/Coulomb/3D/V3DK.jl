export V3DK
"""
Naked 3D Coulomb potential in reciprocal space.
"""
function V3DK(; ϵ = 1, Ω = 1)::Function

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	#Energy unit is eV
	T = qₑ * 1e3 / (ϵ * ϵ₀ * Ω)

	#The unit of input r is Å.
	V = function (k)
		k2 = sum(abs2, k)
		return k2 == 0 ? 0 : T / k2
	end

	return V
end
