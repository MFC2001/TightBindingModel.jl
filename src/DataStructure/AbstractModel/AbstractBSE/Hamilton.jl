
function BSE_NP_Hamilton(vckmap, bandk, bandkq, Kᵈ_F, Kˣ_F, scissor)

	N = length(vckmap)

	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	return BSE_NP_Hamilton!(Htriplet, Hsinglet, vckmap, bandk, bandkq, Kᵈ_F, Kˣ_F, scissor)
end
function BSE_NP_Hamilton!(Htriplet, Hsinglet, vckmap, bandk, bandkq, Kᵈ_F, Kˣ_F, scissor)

	N = length(vckmap)
	tasks = Vector{Task}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 2:N, i in 1:j-1
		n += 1
		tasks[n] = Threads.@spawn begin
			(v′, c′, k′) = vckmap[i]
			(v, c, k) = vckmap[j]

			Kᵈ = Kᵈ_F(v′, c′, k′, v, c, k)
			Kˣ = Kˣ_F(v′, c′, k′, v, c, k)

			Htriplet[i, j] = Kᵈ
			Hsinglet[i, j] = Kᵈ + 2 * Kˣ
		end
	end
	for i in 1:N
		n += 1
		tasks[n] = Threads.@spawn begin
			(v, c, k) = vckmap[i]

			Kᵈ = Kᵈ_F(v, c, k, v, c, k)
			Kˣ = Kˣ_F(v, c, k, v, c, k)

			Δϵ = bandkq[k].values[c] - bandk[k].values[v] + scissor
			Htriplet[i, i] = Δϵ + Kᵈ
			Hsinglet[i, i] = Δϵ + Kᵈ + 2 * Kˣ
		end
	end
	wait.(tasks)

	H_t = Hermitian(Htriplet, :U)
	H_s = Hermitian(Hsinglet, :U)

	return H_t, H_s
end


function BSE_SP_Hamilton(vckmap, bandk, bandkq, Kᵈ_F, Kˣ_F, scissor)

	N = length(vckmap)
	Hbse = Matrix{ComplexF64}(undef, N, N)

	return BSE_SP_Hamilton!(Hbse, vckmap, bandk, bandkq, Kᵈ_F, Kˣ_F, scissor)
end

function BSE_SP_Hamilton!(Hbse, vckmap, bandk, bandkq, Kᵈ_F, Kˣ_F, scissor)

	N = length(vckmap)
	tasks = Vector{Task}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 2:N, i in 1:j-1
		n += 1
		tasks[n] = Threads.@spawn begin
			(v′, c′, k′) = vckmap[i]
			(v, c, k) = vckmap[j]

			Kᵈ = Kᵈ_F(v′, c′, k′, v, c, k)
			Kˣ = Kˣ_F(v′, c′, k′, v, c, k)

			Hbse[i, j] = Kᵈ + Kˣ
		end
	end
	for i in 1:N
		n += 1
		tasks[n] = Threads.@spawn begin
			(v, c, k) = vckmap[i]

			Kᵈ = Kᵈ_F(v, c, k, v, c, k)
			Kˣ = Kˣ_F(v, c, k, v, c, k)

			Δϵ = bandkq[k].values[c] - bandk[k].values[v] + scissor
			Hbse[i, i] = Δϵ + Kᵈ + Kˣ
		end
	end
	wait.(tasks)

	H = Hermitian(Hbse, :U)

	return H
end




# function ()

# Threads.@threads for (i, j) in ittr
# 	(v′, c′, k′) = vck[i]
# 	(v, c, k) = vck[j]

# 	Kᵈ = Kᵈ_F(v′, c′, k′, v, c, k)
# 	Kˣ = Kˣ_F(v′, c′, k′, v, c, k)

# 	if i == j
# 		Δϵ = bandkq[k].values[c] - bandk[k].values[v] + scissor
# 		Htriplet[i, j] = Δϵ + Kᵈ
# 		Hsinglet[i, j] = Δϵ + Kᵈ + 2 * Kˣ
# 	else
# 		Htriplet[i, j] = Kᵈ
# 		Hsinglet[i, j] = Kᵈ + 2 * Kˣ
# 	end
# end
# end



# function BSE_NP_Hamilton!(Htriplet, Hsinglet, ittr, vckindex, bandk, kqindex, Kᵈ_F, Kˣ_F, scissor)

# 	Threads.@threads for (i, j) in ittr
# 		(v₁, c₁, k₁) = vckindex[i]
# 		(v₂, c₂, k₂) = vckindex[j]

# 		Kᵈ = -Kᵈ_F(c₁, kqindex[k₁], v₂, k₂, v₁, k₁, c₂, kqindex[k₂])
# 		Kˣ = Kˣ_F(c₁, kqindex[k₁], v₂, k₂, c₂, kqindex[k₂], v₁, k₁)

# 		if i == j
# 			Δϵ = bandk[kqindex[k₁]].values[c₁] - bandk[k₁].values[v₁] + scissor
# 			Htriplet[i, j] = Δϵ + Kᵈ
# 			Hsinglet[i, j] = Δϵ + Kᵈ + 2 * Kˣ
# 		else
# 			Htriplet[i, j] = Kᵈ
# 			Hsinglet[i, j] = Kᵈ + 2 * Kˣ
# 		end
# 	end

# 	return nothing
# end
