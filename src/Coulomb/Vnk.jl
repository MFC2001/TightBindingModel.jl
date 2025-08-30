
function CreatVnk(kgrid, bandk, U::Function, J::Function)::Function

	Nk = length(kgrid)

	kgrid = reducible_kgrid(kgrid)
	minusmap = kgridmap(kgrid, -)
	addmap = kgridmap(kgrid, +)


	Uk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Jk = Vector{Matrix{ComplexF64}}(undef, Nk)
	Threads.@threads for k in 1:Nk
		Uk[k] = U(kgrid.kdirect[k]) / Nk
		Jk[k] = J(kgrid.kdirect[k]) / Nk
	end

	function Vnk(n₁, k₁, n₂, k₂, n₃, k₃, n₄, k₄)
		ψₙ₁ₖ₁ = bandk[k₁].vectors[:, n₁]
		ψₙ₂ₖ₂ = bandk[k₂].vectors[:, n₂]
		ψₙ₃ₖ₃ = bandk[k₃].vectors[:, n₃]
		ψₙ₄ₖ₄ = bandk[k₄].vectors[:, n₄]
		U₃₂ = Uk[minusmap[k₃, k₂]]
		J₄₂ = Jk[minusmap[k₄, k₂]]
		J₃₄ = Jk[addmap[k₃, k₄]]
		return sum(CartesianIndices(U₃₂)) do I
			(i, j) = Tuple(I)
			return conj(ψₙ₁ₖ₁[i] * ψₙ₂ₖ₂[j]) * (ψₙ₃ₖ₃[j] * ψₙ₄ₖ₄[i] * U₃₂[I] + ψₙ₃ₖ₃[i] * ψₙ₄ₖ₄[j] * J₄₂[I]) +
				   conj(ψₙ₁ₖ₁[i] * ψₙ₂ₖ₂[i]) * ψₙ₃ₖ₃[j] * ψₙ₄ₖ₄[j] * J₃₄[I]
		end
	end

	return Vnk
end

