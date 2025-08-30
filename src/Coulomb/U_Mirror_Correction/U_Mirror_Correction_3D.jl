export U_Mirror_Correction_3D

function U_Mirror_Correction_3D(TB::AbstractTightBindModel, U::HR; αrcut = 4.5, δ = 1e-6, ϵ = 1, head = nothing)

	rcut = maximum(path -> norm(TB.lattice * (path[1:3] + TB.orb_location[path[5]] - TB.orb_location[path[4]])), eachrow(U.path))
	α = αrcut / (rcut * 0.5)

	φR = Gauss_VR(; ϵ, α)

	Ω = abs((TB.lattice[:, 1] × TB.lattice[:, 2]) ⋅ TB.lattice[:, 3])
	φK = Gauss_VK_3D(; ϵ, α, Ω)

	kgrid_max = map(maximum, eachcol(U.path[:, 1:3]))
	kgrid_min = map(minimum, eachcol(U.path[:, 1:3]))
	kgrid = MonkhorstPack((kgrid_max - kgrid_min) .+ 1)
	kgrid = RedKgrid(kgrid)
	Nk = length(kgrid)

	G2_max = -4 * α^2 * log(δ) * 2

	rlattice = reciprocal(TB.lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = abs((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)

	Ggrid = gridindex(Int.(cld.(√G2_max * 1.1, [h₁, h₂, h₃])) * 2 .+ 1)
	Ggrid = filter(G -> sum(abs2, rlattice * G) < G2_max, Ggrid)

	allkG = [k + G for k in kgrid.kdirect, G in Ggrid]

	φkG = map(kG -> φK(rlattice * kG), allkG)


	if isnothing(head)
		qₑ = 1.602176634
		ϵ₀ = 8.854187817
		qsz = (6 * π^2 / (Nk * Ω))^(1 // 3)
		head = 1e3 * qₑ * α * erf(qsz / (2 * α)) / (2 * π^(3 // 2) * ϵ * ϵ₀)
	end

	φ_LR(r) = sum(i -> φkG[i] * cos(2π * (allkG[i] ⋅ r)), eachindex(allkG)) / Nk + head

	correction = Vector{Float64}(undef, size(U.path, 1))
	Threads.@threads for i in axes(U.path, 1)
		path = Vec3(U.path[i, 1], U.path[i, 2], U.path[i, 3])
		i_idx = U.path[i, 4]
		j_idx = U.path[i, 5]
		r_frac = path + TB.orb_location[j_idx] - TB.orb_location[i_idx]
		correction[i] = φ_LR(r_frac) - φR(TB.lattice * r_frac)
	end

	U′ = HR(deepcopy(U.path), U.value - correction; hrsort = 'N', buildhop = 'Y')

	return U′
end

function U_Mirror_Correction_3D_Anisotropy(TB::AbstractTightBindModel, U::HR; αrcut = 4.5, δ = 1e-6, ϵ = I, head = 0)

	error("Can't calculate interaction in real space.")

	rcut = maximum(path -> norm(TB.lattice * (path[1:3] + TB.orb_location[path[5]] - TB.orb_location[path[4]])), eachrow(U.path))
	α = αrcut / (rcut * 0.5)


	Ω = abs((TB.lattice[:, 1] × TB.lattice[:, 2]) ⋅ TB.lattice[:, 3])
	φK = Gauss_VK_3D_Anisotropy(; ϵ, α, Ω)

	Tϵ = eigen(ϵ)

	#1e-19
	qₑ = 1.602176634
	#1e-12
	ϵ₀ = 8.854187817
	#Energy unit is eV
	T = qₑ * 1e3 / (prod(Tϵ.values) * ϵ₀ * (2π)^3)


	A(sθ, cθ, sφ, cφ) = sθ^2 * cφ^2 / Tϵ.values[1] + sθ^2 * sφ^2 / Tϵ.values[2] + cθ^2 / Tϵ.values[3]

	Nθ = 200
	dk = 0.01

	dkθφ = dk * π^2 / Nθ^2
	θ = range(0, π, Nθ)
	Nφ = 2Nθ
	φ = range(0, 2π, Nφ)
	sθ = sin.(θ)
	cθ = cos.(θ)
	sφ = sin.(φ)
	cφ = cos.(φ)
	Aθφ = [A(sθ[θi], cθ[θi], sφ[φi], cφ[φi]) for θi in 1:Nθ, φi in 1:Nφ]

	k_max = √(-4 * α^2 * log(δ) / minimum(Aθφ))
	Nk = Int(cld(k_max, dk))
	k = range(0, step = dk, length = Nk)

	T = T * dkθφ
	sqrtϵ = sqrt.(Tϵ.values)
	return φR = function (r)
		r = norm(Tϵ.vectors * r ./ sqrtϵ)
		result = 0
		for ki in 1:Nk
			k2α = k[ki]^2 / (4 * α)
			for θi in 1:Nθ, φi in 1:Nφ
				result += sθ[θi] * cis(k[ki] * r * cθ[θi]) * exp(-Aθφ[θi, φi] * k2α)
			end
		end
		return result * T
	end


	kgrid_max = map(maximum, eachcol(U.path[:, 1:3]))
	kgrid_min = map(minimum, eachcol(U.path[:, 1:3]))
	kgrid_size = (kgrid_max - kgrid_min) .+ 1

	# Ggrid
	G2_max = -4 * α^2 * log(δ) * 2

	rlattice = reciprocal(TB.lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = abs((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)

	Ggrid = gridindex(Int.(cld.(√G2_max * 1.1, [h₁, h₂, h₃])) * 2 .+ 1)
	Ggrid = filter(G -> sum(abs2, rlattice * G) < G2_max, Ggrid)


	kgrid = MonkhorstPack(kgrid_size)
	kgrid = RedKgrid(kgrid)
	Nk = length(kgrid)

	allkG = [k + G for k in kgrid.kdirect, G in Ggrid]
	φkG = map(kG -> φK(rlattice * kG), allkG)
	return φ_LR(r) = sum(i -> φkG[i] * cis(2π * (allkG[i] ⋅ r)), eachindex(allkG)) / Nk + head



	correction = map(eachrow(U.path)) do path
		r_frac = path[1:3] + TB.orb_location[path[5]] - TB.orb_location[path[4]]
		Threads.@spawn φ_LR(r_frac) - φR(TB.lattice * r_frac)
	end
	correction = fetch.(correction)

	U′ = HR(U.path, U.value - correction; hrsort = 'N', buildhop = 'Y')

	return U′
end
