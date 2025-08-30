export U_Mirror_Correction_2D

function U_Mirror_Correction_2D(TB::AbstractTightBindModel, U::HR; αrcut = 4.5, δ = 1e-6, ϵ = 1)

	rcut = maximum(path -> norm(TB.lattice * (path[1:3] + TB.orb_location[path[5]] - TB.orb_location[path[4]])), eachrow(U.path))
	α = αrcut / (rcut * 0.5)

	φR = Gauss_VR(; ϵ, α)

	#kgrid
	kgrid_max = map(maximum, eachcol(U.path[:, 1:3]))
	kgrid_min = map(minimum, eachcol(U.path[:, 1:3]))
	kgrid = MonkhorstPack((kgrid_max - kgrid_min) .+ 1)
	kgrid = RedKgrid(kgrid)
	Nk = length(kgrid)

	#Ggrid
	G2_max = -4 * α^2 * log(δ) * 2

	rlattice = reciprocal(TB.lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = abs((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)
	Ggrid = Int.(cld.(√G2_max * 1.1, [h₁, h₂, h₃])) * 2 .+ 1

	a₁ = TB.lattice[:, 1]
	a₂ = TB.lattice[:, 2]
	a₃ = TB.lattice[:, 3]

	if TB.period[1] == "np"
		Ggrid[1] = 1
		S = norm(a₂ × a₃)
		zindex = 1
	elseif TB.period[2] == "np"
		Ggrid[2] = 1
		S = norm(a₃ × a₁)
		zindex = 2
	elseif TB.period[3] == "np"
		Ggrid[3] = 1
		S = norm(a₁ × a₂)
		zindex = 3
	else
		error("Wrong TB.period.")
	end
	φK = Gauss_VK_2D(; ϵ, α, S)


	Ggrid = gridindex(Ggrid)
	Ggrid = filter(G -> sum(abs2, rlattice * G) < G2_max, Ggrid)

	allkG = reshape([k + G for k in kgrid.kdirect, G in Ggrid], :)


	# ϕ₀
	qₑ = 1.602176634
	ϵ₀ = 8.854187817
	qsz = √(4π / (Nk * S))
	T = 4π * ϵ * ϵ₀
	qα = qsz / (2 * α)
	φ₀_z0 = 1e3 * qₑ * (qsz * erfc(qα) - 2 * α * exp(-(qsz / 2 * α)^2) / √π + 2 * α / √π) / T

	φ₀_z = function (z)
		z = abs(z)
		if iszero(z)
			return φ₀_z0
		else
			αz = α * z
			eqz = exp(qsz * z)
			return 1e3 * qₑ * (erfc(qα + αz) * eqz - erfc(qα - αz) / eqz + 2 * erf(αz)) / (2 * T * z)
		end
	end



	norb = length(TB.orb_location)
	φkG = Array{Float64}(undef, length(allkG), norb, norb)
	φ₀ = Matrix{Float64}(undef, norb, norb)

	allkG_car = map(kG -> rlattice * kG, allkG)
	for i in 1:norb, j in 1:i
		dorb = TB.lattice * (TB.orb_location[j] - TB.orb_location[i])
		z = dorb[zindex]
		φkG[:, i, j] = map(kG -> φK(kG, z), allkG_car)
		φkG[:, j, i] = φkG[:, i, j]
		φ₀[i, j] = φ₀_z(z)
		φ₀[j, i] = φ₀[i, j]
	end

	I = findfirst(iszero, allkG)
	φkG[I, :, :] .= 0


	φ_LR(r, i, j) = sum(kGi -> φkG[kGi, i, j] * cos(2π * (allkG[kGi] ⋅ r)), eachindex(allkG)) / Nk + φ₀[i, j]
	# φ_LR(r, i, j) = sum(kGi -> φkG[kGi, i, j] * cis(2π * (allkG[kGi] ⋅ r)), eachindex(allkG)) / Nk + head


	correction = Vector{Float64}(undef, size(U.path, 1))
	Threads.@threads for i in axes(U.path, 1)
		path = Vec3(U.path[i, 1], U.path[i, 2], U.path[i, 3])
		i_idx = U.path[i, 4]
		j_idx = U.path[i, 5]
		r_frac = path + TB.orb_location[j_idx] - TB.orb_location[i_idx]
		correction[i] = φ_LR(r_frac, i_idx, j_idx) - φR(TB.lattice * r_frac)
	end

	U′ = HR(U.path, U.value - correction; hrsort = 'N', buildhop = 'Y')

	return U′
end
