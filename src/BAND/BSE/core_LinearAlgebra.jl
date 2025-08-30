function BSE_NP(qpoints, bse::AbstractBSE, ::Val{true})

	Nq = length(qpoints)
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)

		(bandkq, Kᵈ, Kˣ) = Kernal!(bse, q)
		(H_t, H_s) = BSE_NP_Hamilton!(Htriplet, Hsinglet, bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)

		BSEband_t[qi] = eigen!(H_t)
		BSEband_s[qi] = eigen!(H_s)
	end

	return BSEband_t, BSEband_s
end
function BSE_NP(qpoints, bse::AbstractBSE, ::Val{false})

	Nq = length(qpoints)
	N = length(bse.vckmap)
	BSEband_t = Matrix{Float64}(undef, N, Nq)
	BSEband_s = Matrix{Float64}(undef, N, Nq)


	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)

		(bandkq, Kᵈ, Kˣ) = Kernal!(bse, q)
		(H_t, H_s) = BSE_NP_Hamilton!(Htriplet, Hsinglet, bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)

		BSEband_t[:, qi] = eigvals!(H_t)
		BSEband_s[:, qi] = eigvals!(H_s)
	end

	return BSEband_t, BSEband_s
end
function BSE_SP(qpoints, bse::AbstractBSE, ::Val{true})

	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	H_bse = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)

		(bandkq, Kᵈ, Kˣ) = Kernal!(bse, q)
		Hbse = BSE_SP_Hamilton!(H_bse, bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)

		BSEband[qi] = eigen!(Hbse)
	end

	return BSEband
end
function BSE_SP(qpoints, bse::AbstractBSE, ::Val{false})

	Nq = length(qpoints)
	N = length(bse.vckmap)
	BSEband = Matrix{Float64}(undef, N, Nq)

	N = length(bse.vckmap)
	H_bse = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)

		(bandkq, Kᵈ, Kˣ) = Kernal!(bse, q)
		Hbse = BSE_SP_Hamilton!(H_bse, bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)

		BSEband[:, qi] = eigvals!(Hbse)
	end

	return BSEband
end
