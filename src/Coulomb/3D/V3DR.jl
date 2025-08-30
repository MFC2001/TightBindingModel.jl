export V3DR
"""
	V3DR(; ϵ = 1, V₀ = 0, r₀ = 1e-3)::Function
	
Naked 3D Coulomb potential in real space.
"""
function V3DR(; ϵ = 1)::Function

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	#Energy unit is eV
	T = qₑ * 1e3 / (4 * π * ϵ * ϵ₀)

	#The unit of input r is Å.
	V = function (r)
		r = norm(r)
		return r == 0 ? 0 : T / r
	end

	return V
end
