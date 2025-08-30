export Gauss_VK_2D, Gauss_VK_3D, Gauss_VK_3D_Anisotropy


function Gauss_VK_3D(; ϵ = 1, α = 1, Ω = 1)::Function

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	#Energy unit is eV
	T = qₑ * 1e3 / (ϵ * ϵ₀ * Ω)

	α²4 = 4 * α^2

	#The unit of input r is Å.
	V = function (k)
		k2 = sum(abs2, k)
		return k2 == 0 ? 0 : T * exp(-k2 / α²4) / k2
	end

	return V
end

function Gauss_VK_2D(; ϵ = 1, α = 1, S = 1)::Function

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	#Energy unit is eV
	T = qₑ * 1e3 / (4 * ϵ * ϵ₀ * S)

	α2 = 2 * α

	#The unit of input r is Å.
	V = function (k, z)
		k = norm(k)
		z = abs(z)

		# ekz = exp(k * z)
		# k2α = k / (α2)
		# αz = α * z

		# return k == 0 ? 0 : T * (erfc(k2α + αz) * ekz + erfc(k2α - αz) / ekz) / k


		if k == 0
			αz = α * z
			Vk = T * z * (erfc(αz) - erfc(-αz))
		else
			ekz = exp(k * z)
			k2α = k / (α2)
			αz = α * z

			Vk = T * (erfc(k2α + αz) * ekz + erfc(k2α - αz) / ekz) / k
		end

		return Vk
	end

	return V
end

function Gauss_VK_3D_Anisotropy(; ϵ = I, α = 1, Ω = 1)::Function

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817

	#Energy unit is eV
	T = qₑ * 1e3 / (ϵ₀ * Ω)

	α²4 = 4 * α^2

	#The unit of input r is Å.
	V = function (k)
		k2 = sum(abs2, k)
		return k2 == 0 ? 0 : T * exp(-k2 / α²4) / (transpose(k) * ϵ * k)
	end

	return V
end
