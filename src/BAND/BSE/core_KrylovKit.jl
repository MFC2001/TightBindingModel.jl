function BSE_NP(qpoints, bse::AbstractBSE, ::Val{true})

	Nq = length(qpoints)
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)

		bandkq, Kᵈ, Kˣ = Kernal!(bse, q)
		H_t, H_s = BSE_NP_Hamilton!(Htriplet, Hsinglet, bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)

		BSEband_t[qi] = _eigsolve_Hmat(H_t)
		BSEband_s[qi] = _eigsolve_Hmat(H_s)
	end

	return BSEband_t, BSEband_s
end
function BSE_NP(qpoints, bse::AbstractBSE, ::Val{false})

	(BSEband_t, BSEband_s) = BSE_NP(qpoints, bse, Val(true))

	BSEband_t = _eigen2vals(BSEband_t)
	BSEband_s = _eigen2vals(BSEband_s)

	return BSEband_t, BSEband_s
end

function BSE_SP(qpoints, bse::AbstractBSE, ::Val{true})

	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	H_bse = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)

		bandkq, Kᵈ, Kˣ = Kernal!(bse, q)
		Hbse = BSE_SP_Hamilton!(H_bse, bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)

		BSEband[qi] = _eigsolve_Hmat(Hbse)
	end

	return BSEband
end
function BSE_SP(qpoints, bse::AbstractBSE, ::Val{false})

	BSEband = BSE_SP(qpoints, bse, Val(true))
	BSEband = _eigen2vals(BSEband)

	return BSEband
end
