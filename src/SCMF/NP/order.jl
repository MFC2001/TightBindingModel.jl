
function MForder_NP(K, ρ)
	Threads.@threads for i in axes(K, 1)
		Δ[i] = sum(I -> K[i, I] * ρ[I], CartesianIndices(ρ))
	end
	return Δ
end

function SCMF_NP_kernal(ViR, Δindex, R, norb)

	nΔ = length(Δindex)
	NR = length(R)

	addmap = gridmap(R, +)
	ΓR = findfirst(iszero, R)


	K = Array{ComplexF64}(undef, nΔ, norb, norb, NR)
	Threads.@threads for R′ in 1:NR
		for n in 1:nΔ, i in 1:norb, j in 1:norb
			(α, β, R) = Δindex[n]
			K[n, i, j, R′] = sum(1:NR) do R₀′
				R₁′ = addmap[R₀′, R′]
				# return ViR(i, R₀′, α, ΓR, β, R, j, R₁′;spin=true) - ViR(i, R₀′, α, ΓR, j, R₁′, β, R;spin=true)
				return ViR(i, R₀′, α, ΓR, β, R, j, R₁′; spin = true) - ViR(i, R₀′, α, ΓR, j, R₁′, β, R; spin = true) + ViR(j, R₁′, α, ΓR, β, R, i, R₀′; spin = false)
			end
		end
	end

	return K
end
