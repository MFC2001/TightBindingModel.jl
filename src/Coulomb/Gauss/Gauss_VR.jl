export Gauss_VR

function Gauss_VR(; ϵ = 1, α = 1)::Function

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	#Energy unit is eV
	T = qₑ * 1e3 / (4 * π * ϵ * ϵ₀)

	V₀ = T * 2 * α / √π

	#The unit of input r is Å.
	V = function (r)
		r = norm(r)
		return r == 0 ? V₀ : T * erf(α * r) / r
	end

	return V
end

# function Gauss_VR_Anisotropy


# 	Vk = Gauss_VK_3D_Anisotropy(; ϵ = I, α = 1, Ω = 1)
# end