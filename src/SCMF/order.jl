


function MForder!(Δ, Δittr, ijkittr, Vik, ρ)

	Threads.@threads for (α, β, k) in Δittr
		Δ[k][α, β] = sum(ijkittr) do (i, j, k′)
			(Vik(α, k, i, k′, j, k′, β, k) - Vik(α, k, i, k′, β, k, j, k′)) * ρ[k′][i, j]
		end
		Δ[k][β, α] = conj(Δ[k][α, β])
	end

	return Δ
end



function generate_independent_order(lattsymops)
	norb = length(lattsymops)
	Δ = Array{ComplexF64}(undef,norb,norb,Nk)

	


end