using SpecialFunctions
using Struve
function V2DK_R(V₀::Real; ϵ, LC, ϵ₁ = 1, ϵ₂ = 1, r₀ = 1)

	if ϵ == 1
		Error("When ϵ=1, it's equal to the bare Coulomb potential in 3D.")
	end

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	ρ₀ = LC * (ϵ - 1) / 2
	T = qₑ * 1e3 / (8 * ϵ₀ * ρ₀)

	ϵd = (ϵ₁ + ϵ₂) / 2

	V = function (r)
		r = norm(r)
		x = ϵd * r / ρ₀
		return r < r₀ ? V₀ : T * (struveh(0, x) - bessely0(x))
	end

	return V
end
