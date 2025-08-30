
function BSE_NP_A(
	bse::BSE,
	qindex::Integer;
	scissor = 0,
	vector::Bool = false,
	arpack::Bool = false,
	NE = 10,
	Emin = 0,
)

	nv = length(bse.vindex)
	nc = length(bse.cindex)
	Nk = length(bse.kgrid)

	N = nv * nc * Nk
	vckindex = Vector{Tuple{Int, Int, Int}}(undef, N)
	n = 0
	for k in 1:Nk, c in bse.cindex, v in bse.vindex
		n += 1
		vckindex[n] = (v, c, k)
	end

	ittr = Vector{Tuple{Int, Int}}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 1:N, i in 1:j
		n += 1
		ittr[n] = (i, j)
	end

	kgrid = reducible_kgrid(bse.kgrid)
	addmap = kgridmap(kgrid, +)
	minusmap = kgridmap(kgrid, -)

	bandk = BAND(kgrid, bse.TB; vector = true)


	U_screened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	J_screened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	A_screened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	U_unscreened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	J_unscreened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	A_unscreened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	Threads.@threads for k in 1:Nk
		U_screened_k[k] = bse.U_screened(kgrid.kdirect[k]) / Nk
		J_screened_k[k] = bse.J_screened(kgrid.kdirect[k]) / Nk
		A_screened_k[k] = bse.A_screened(kgrid.kdirect[k]) / Nk
		U_unscreened_k[k] = bse.U_unscreened(kgrid.kdirect[k]) / Nk
		J_unscreened_k[k] = bse.J_unscreened(kgrid.kdirect[k]) / Nk
		A_unscreened_k[k] = bse.A_unscreened(kgrid.kdirect[k]) / Nk
	end


	if qindex == 0
		kqindex = collect(1:Nk)
	else
		addmap = kgridmap(kgrid, +)
		kqindex = addmap[:, q]
	end


	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)


	bandkq = bandk[kqindex]
	J_screened_kq = J_screened_k[kqindex]
	A_screened_kq = A_screened_k[kqindex]
	U_unscreened_kq = U_unscreened_k[kqindex]
	J_unscreened_kq = J_unscreened_k[kqindex]
	A_unscreened_kq = A_unscreened_k[kqindex]


	Kᵈ = CreatKernal_A(bandkq, bandk, bandk, bandkq, U_screened_k, J_screened_kq, J_screened_kq,
		A_screened_kq, A_screened_k, A_screened_k, A_screened_kq, addmap, minusmap)
	Kˣ = CreatKernal_A(bandkq, bandk, bandkq, bandk, U_unscreened_kq, J_unscreened_k, J_unscreened_kq,
		A_unscreened_kq, A_unscreened_k, A_unscreened_kq, A_unscreened_k, addmap, minusmap)


	BSE_NP_Hamilton!(Htriplet, Hsinglet, ittr, vckindex, bandk, bandkq, Kᵈ, Kˣ, scissor)
	H_t = Hermitian(Htriplet, :U)
	H_s = Hermitian(Hsinglet, :U)


	if arpack
		H_t = eigs(H_t; nev = NE, sigma = Emin, ritzvec = vector)
		H_s = eigs(H_s; nev = NE, sigma = Emin, ritzvec = vector)
		E_t = real(H_t[1])
		I_t = sortperm(E_t)
		E_s = real(H_s[1])
		I_s = sortperm(E_s)
		if vector
			BSEband_t = Eigen(E_t[I_t], H_t[2][:, I_t])
			BSEband_s = Eigen(E_s[I_s], H_s[2][:, I_s])
		else
			BSEband_t = E_t[I_t]
			BSEband_s = E_s[I_s]
		end
	else
		if vector
			BSEband_t = eigen!(H_t)
			BSEband_s = eigen!(H_s)
		else
			BSEband_t = eigvals!(H_t)
			BSEband_s = eigvals!(H_s)
		end
	end



	if vector
		wave = Vector{Function}(undef, 4)

		wave[1] = BSE_NP_Q_wave_fixe(BSEband_t, bandk, kgrid, vckindex, kqindex)
		wave[2] = BSE_NP_Q_wave_fixh(BSEband_t, bandk, kgrid, vckindex, kqindex)
		wave[3] = BSE_NP_Q_wave_fixe(BSEband_s, bandk, kgrid, vckindex, kqindex)
		wave[4] = BSE_NP_Q_wave_fixh(BSEband_s, bandk, kgrid, vckindex, kqindex)

		return BSEband_t, BSEband_s, wave
	else
		return BSEband_t, BSEband_s
	end

end

function BSE_NP_A(
	bse::BSE,
	kline::Kline;
	scissor = 0,
	vector::Bool = false,
	arpack::Bool = false,
	NE = 10,
	Emin = 0,
)

	nv = length(bse.vindex)
	nc = length(bse.cindex)
	Nk = length(bse.kgrid)
	Nq = length(kline)


	N = nv * nc * Nk
	vckindex = Vector{Tuple{Int, Int, Int}}(undef, N)
	n = 0
	for k in 1:Nk, c in bse.cindex, v in bse.vindex
		n += 1
		vckindex[n] = (v, c, k)
	end

	ittr = Vector{Tuple{Int, Int}}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 1:N, i in 1:j
		n += 1
		ittr[n] = (i, j)
	end


	kgrid = reducible_kgrid(bse.kgrid)
	addmap = kgridmap(kgrid, +)
	minusmap = kgridmap(kgrid, -)

	bandk = BAND(kgrid, bse.TB; vector = true)


	U_screened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	A_screened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	J_unscreened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	A_unscreened_k = Vector{Matrix{ComplexF64}}(undef, Nk)
	Threads.@threads for k in 1:Nk
		U_screened_k[k] = bse.U_screened(kgrid.kdirect[k]) / Nk
		A_screened_k[k] = bse.A_screened(kgrid.kdirect[k]) / Nk
		J_unscreened_k[k] = bse.J_unscreened(kgrid.kdirect[k]) / Nk
		A_unscreened_k[k] = bse.A_unscreened(kgrid.kdirect[k]) / Nk
	end



	J_screened_kq = Vector{Matrix{ComplexF64}}(undef, Nk)
	A_screened_kq = Vector{Matrix{ComplexF64}}(undef, Nk)
	U_unscreened_kq = Vector{Matrix{ComplexF64}}(undef, Nk)
	J_unscreened_kq = Vector{Matrix{ComplexF64}}(undef, Nk)
	A_unscreened_kq = Vector{Matrix{ComplexF64}}(undef, Nk)

	if vector
		BSEband_t = Vector{Eigen}(undef, Nq)
		BSEband_s = Vector{Eigen}(undef, Nq)
	else
		if arpack && NE < N
			BSEband_t = Matrix{Float64}(undef, NE, Nq)
			BSEband_s = Matrix{Float64}(undef, NE, Nq)
		else
			BSEband_t = Matrix{Float64}(undef, N, Nq)
			BSEband_s = Matrix{Float64}(undef, N, Nq)
		end
	end

	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)


	for (qi, q) in enumerate(kline.kdirect)

		bandkq = BAND(kgrid.kdirect .+ [q], bse.TB; vector = true)

		Threads.@threads for k in 1:Nk
			J_screened_kq[k] = bse.J_screened(kgrid.kdirect[k] + q) / Nk
			A_screened_kq[k] = bse.A_screened(kgrid.kdirect[k] + q) / Nk
			U_unscreened_kq[k] = bse.U_unscreened(kgrid.kdirect[k] + q) / Nk
			J_unscreened_kq[k] = bse.J_unscreened(kgrid.kdirect[k] + q) / Nk
			A_unscreened_kq[k] = bse.A_unscreened(kgrid.kdirect[k] + q) / Nk
		end


		Kᵈ = CreatKernal_A(bandkq, bandk, bandk, bandkq, U_screened_k, J_screened_kq, J_screened_kq,
			A_screened_kq, A_screened_k, A_screened_k, A_screened_kq, addmap, minusmap)
		Kˣ = CreatKernal_A(bandkq, bandk, bandkq, bandk, U_unscreened_kq, J_unscreened_k, J_unscreened_kq,
			A_unscreened_kq, A_unscreened_k, A_unscreened_kq, A_unscreened_k, addmap, minusmap)

		BSE_NP_Hamilton!(Htriplet, Hsinglet, ittr, vckindex, bandk, bandkq, Kᵈ, Kˣ, scissor)
		H_t = Hermitian(Htriplet, :U)
		H_s = Hermitian(Hsinglet, :U)

		if arpack
			H_t = eigs(H_t; nev = NE, sigma = Emin, ritzvec = vector)
			H_s = eigs(H_s; nev = NE, sigma = Emin, ritzvec = vector)
			E_t = real(H_t[1])
			I_t = sortperm(E_t)
			E_s = real(H_s[1])
			I_s = sortperm(E_s)
			if vector
				BSEband_t[qi] = Eigen(E_t[I_t], H_t[2][:, I_t])
				BSEband_s[qi] = Eigen(E_s[I_s], H_s[2][:, I_s])
			else
				BSEband_t[:, qi] = E_t[I_t]
				BSEband_s[:, qi] = E_s[I_s]
			end
		else
			if vector
				BSEband_t[qi] = eigen!(H_t)
				BSEband_s[qi] = eigen!(H_s)
			else
				BSEband_t[:, qi] = eigvals!(H_t)
				BSEband_s[:, qi] = eigvals!(H_s)
			end
		end

	end


	if vector
		wave = Vector{Function}(undef, 4)

		wave[1] = BSE_NP_wave_fixe(BSEband_t, bandk, kgrid, vckindex, addmap)
		wave[2] = BSE_NP_wave_fixh(BSEband_t, bandk, kgrid, vckindex, addmap)
		wave[3] = BSE_NP_wave_fixe(BSEband_s, bandk, kgrid, vckindex, addmap)
		wave[4] = BSE_NP_wave_fixh(BSEband_s, bandk, kgrid, vckindex, addmap)

		return BSEband_t, BSEband_s, wave
	else
		return BSEband_t, BSEband_s
	end

end