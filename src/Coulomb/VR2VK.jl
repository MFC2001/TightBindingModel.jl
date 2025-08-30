export VR2VK
"""
	LatticeCoulomb(VR::HR; Γ = 1)::Function

	Vᵢⱼ(k) = ∑_{R} U_{i0jR} exp(ik⋅R)

	Return a function that can be used by fun(k::AbstractVector{<:Real}) to get Vmatrix at this kpoint.
"""
function VR2VK(VR::HR)

	norb = numorb(VR)

	# input direct reciprocical coordinate.
	function Vk!(A, k)
		size(A) == (norb, norb) || error("Buffer size mismatch.")
		A .= map(CartesianIndices(VR.Nhop)) do I
			iszero(VR.Nhop[I]) ? 0 : sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), VR.hop[I])
		end
		return A
	end

	function Vk(k)
		return map(CartesianIndices(VR.Nhop)) do I
			iszero(VR.Nhop[I]) ? 0 : sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), VR.hop[I])
		end
	end

	return (f = Vk, f! = Vk!)
end

"""
	LatticeCoulomb(VR::HR, orblocat::AbstractVector; Γ = 1)::Function

	Vᵢⱼ(k) = ∑_{R} U_{i0jR} exp(ik⋅(R+τⱼ-τᵢ))

	Return a function that can be used by fun(k::AbstractVector{<:Real}) to get Vmatrix at this kpoint.
"""
function VR2VK(VR::HR, orblocat::AbstractVector{<:ReducedCoordinates{<:Real}})

	norb = numorb(VR)

	function Vk!(A, k)
		size(A) == (norb, norb) || error("Buffer size mismatch.")
		A .= map(CartesianIndices(VR.Nhop)) do I
			(i, j) = Tuple(I)
			iszero(VR.Nhop[I]) ? 0 : cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), VR.hop[I])
		end
	end

	function Vk(k)
		return map(CartesianIndices(VR.Nhop)) do I
			(i, j) = Tuple(I)
			iszero(VR.Nhop[I]) ? 0 : cis(2π * (k ⋅ (orblocat[j] - orblocat[i]))) * sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), VR.hop[I])
		end
	end

	return (f = Vk, f! = Vk!)
end
