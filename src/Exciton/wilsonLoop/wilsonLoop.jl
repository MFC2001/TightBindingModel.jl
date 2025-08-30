export BSE_NP_WilsonLoop

function BSE_NP_Wilsonloop(bse::BSE;
	bandindex_t::Union{Nothing, AbstractVector{<:Integer}, Integer} = nothing,
	bandindex_s::Union{Nothing, AbstractVector{<:Integer}, Integer} = nothing,
	outfolder = "./",
	arpack::Bool = false,
	NE = 10,
	Emin = 0,
)

	qgrid = bse.kgrid
	#basis is twice as dense
	bse.kgrid = MonkhorstPack(qgrid.kgrid_size .* 2; kshift = [1 // 2, 1 // 2, 1 // 2])
	nv = length(bse.v)
	nc = length(bse.c)
	Nk = length(bse.kgrid)
	N = nv * nc * Nk
	bse.vckindex = Vector{Tuple{Int, Int, Int}}(undef, N)
	n = 0
	for k in 1:Nk, c in 1:nc, v in 1:nv
		n += 1
		bse.vckindex[n] = (v, c, k)
	end

	qgrid = MonkhorstPack(qgrid.kgrid_size)
	redqgrid = reducible_kgrid(qgrid)
	(BSEband_t, BSEband_s) = BSE_NP(bse, redqgrid.kdirect; vector = true, arpack, NE, Emin)


	if !isnothing(bandindex_t)
		if bandindex_t isa Integer
			wcc_qx, wcc_qy = _WilsonLoop_BSE_NP_oneband(qgrid, bse, BSEband_t, bandindex_t)
		end
		if length(bandindex_t) == 1
			wcc_qx, wcc_qy = _WilsonLoop_BSE_NP_oneband(qgrid, bse, BSEband_t, bandindex_t[1])
		end
		wcc_qx, wcc_qy = _WilsonLoop_BSE_NP_multiband(qgrid, bse, BSEband_t, bandindex_t)
		WriteBAND(wcc_qx, joinpath(outfolder, "wcc_qx_triplet.dat"))
		WriteBAND(wcc_qy, joinpath(outfolder, "wcc_qy_triplet.dat"))
	end
	if !isnothing(bandindex_s)
		if bandindex_s isa Integer
			wcc_qx, wcc_qy = _WilsonLoop_BSE_NP_oneband(qgrid, bse, BSEband_s, bandindex_s)
		end
		if length(bandindex_s) == 1
			wcc_qx, wcc_qy = _WilsonLoop_BSE_NP_oneband(qgrid, bse, BSEband_s, bandindex_s[1])
		end
		wcc_qx, wcc_qy = _WilsonLoop_BSE_NP_multiband(qgrid, bse, BSEband_s, bandindex_s)
		WriteBAND(wcc_qx, joinpath(outfolder, "wcc_qx_singlet.dat"))
		WriteBAND(wcc_qy, joinpath(outfolder, "wcc_qy_singlet.dat"))
	end

	return nothing
end

function _WilsonLoop_BSE_NP_oneband(qgrid, bse, BSEband, bandindex)

	redqgrid = reducible_kgrid(qgrid)

	Nqx = qgrid.kgrid_size[1]
	Nqy = qgrid.kgrid_size[2]

	start = -floor.(Int, ([Nqx, Nqy] .- 1) .// 2)
	stop = ceil.(Int, ([Nqx, Nqy] .- 1) .// 2)


	#W(qy)
	qy = start[2]-1:stop[2]
	nnkpts = zeros(Int, 8, Nqx)
	θ_qy = Vector{Float64}(undef, length(qy))
	for qyi in eachindex(qy)
		aimqdirects = [Vec3([qx // Nqx, qy[qyi] // Nqy, 0]) for qx in start[1]-1:stop[1]]
		aimq_qindex = Vector{Int}(undef, Nqx + 1)
		for (qi, qv) in enumerate(aimqdirects)
			aimq_qindex[qi] = findfirst(q -> all(isinteger, q - qv), redqgrid.kdirect)
		end
		G = Int.(aimqdirects[1] - redqgrid.kdirect[aimq_qindex[1]])
		nnkpts[5, 1] = aimq_qindex[1]
		nnkpts[6:8, 1] = G
		for i in 2:Nqx
			G = Int.(aimqdirects[i] - redqgrid.kdirect[aimq_qindex[i]])
			nnkpts[1, i-1] = aimq_qindex[i]
			nnkpts[2:4, i-1] = G
			nnkpts[5, i] = aimq_qindex[i]
			nnkpts[6:8, i] = G
		end
		G = Int.(aimqdirects[end] - redqgrid.kdirect[aimq_qindex[end]])
		nnkpts[1, end] = aimq_qindex[end]
		nnkpts[2:4, end] = G

		M = _Wilsonloop_mmn_BSE_NP(qgrid, bse, BSEband, bandindex, nnkpts)

		W = prod(M)
		θ_qy[qyi] = angle(W)
	end
	θ_qy = smooth_AngleSingular(θ_qy, π)
	wcc_qy = θ_qy ./ 2π

	#W(qx)
	qx = start[1]-1:stop[1]
	nnkpts = zeros(Int, 8, Nqy)
	θ_qx = Vector{Float64}(undef, length(qx))
	for qxi in eachindex(qx)
		aimqdirects = [Vec3([qx[qxi] // Nqx, ky // Nqy, 0]) for ky in start[2]-1:stop[2]]
		aimq_qindex = Vector{Int}(undef, Nqy + 1)
		for (qi, qv) in enumerate(aimqdirects)
			aimq_qindex[qi] = findfirst(q -> all(isinteger, q - qv), redqgrid.kdirect)
		end
		G = Int.(aimqdirects[1] - redqgrid.kdirect[aimq_qindex[1]])
		nnkpts[5, 1] = aimq_qindex[1]
		nnkpts[6:8, 1] = G
		for i in 2:Nqy
			G = Int.(aimqdirects[i] - redqgrid.kdirect[aimq_qindex[i]])
			nnkpts[1, i-1] = aimq_qindex[i]
			nnkpts[2:4, i-1] = G
			nnkpts[5, i] = aimq_qindex[i]
			nnkpts[6:8, i] = G
		end
		G = Int.(aimqdirects[end] - redqgrid.kdirect[aimq_qindex[end]])
		nnkpts[1, end] = aimq_qindex[end]
		nnkpts[2:4, end] = G

		M = _Wilsonloop_mmn_BSE_NP(qgrid, bse, BSEband, bandindex, nnkpts)

		W = prod(M)
		θ_qx[qxi] = angle(W)
	end
	θ_qx = smooth_AngleSingular(θ_qx, π)
	wcc_qx = θ_qx ./ 2π

	return wcc_qx, wcc_qy
end
function _WilsonLoop_BSE_NP_multiband(qgrid, bse, BSEband, bandindex)

	redqgrid = reducible_kgrid(qgrid)
	nband = length(bandindex)

	Nqx = qgrid.kgrid_size[1]
	Nqy = qgrid.kgrid_size[2]

	start = -floor.(Int, ([Nqx, Nqy] .- 1) .// 2)
	stop = ceil.(Int, ([Nqx, Nqy] .- 1) .// 2)


	#W(qy)
	qy = start[2]-1:stop[2]
	nnkpts = zeros(Int, 8, Nqx)
	θ_qy = Matrix{Float64}(undef, nband, length(qy))
	for qyi in eachindex(qy)
		aimqdirects = [Vec3([qx // Nqx, qy[qyi] // Nqy, 0]) for qx in start[1]-1:stop[1]]
		aimq_qindex = Vector{Int}(undef, Nqx + 1)
		for (qi, qv) in enumerate(aimqdirects)
			aimq_qindex[qi] = findfirst(q -> all(isinteger, q - qv), redqgrid.kdirect)
		end
		G = Int.(aimqdirects[1] - redqgrid.kdirect[aimq_qindex[1]])
		nnkpts[5, 1] = aimq_qindex[1]
		nnkpts[6:8, 1] = G
		for i in 2:Nqx
			G = Int.(aimqdirects[i] - redqgrid.kdirect[aimq_qindex[i]])
			nnkpts[1, i-1] = aimq_qindex[i]
			nnkpts[2:4, i-1] = G
			nnkpts[5, i] = aimq_qindex[i]
			nnkpts[6:8, i] = G
		end
		G = Int.(aimqdirects[end] - redqgrid.kdirect[aimq_qindex[end]])
		nnkpts[1, end] = aimq_qindex[end]
		nnkpts[2:4, end] = G

		M = _Wilsonloop_mmn_BSE_NP(qgrid, bse, BSEband, bandindex, nnkpts)

		W = 1
		for i in 1:Nqx
			F = svd(M[:, :, i])
			W = (F.U * F.Vt) * W
		end

		θ_qy[:, qyi] = eigvals!(Hermitian(W))
	end
	wcc_qy = θ_qy ./ 2π

	#W(qx)
	qx = start[1]-1:stop[1]
	nnkpts = zeros(Int, 8, Nqy)
	θ_qx = Matrix{Float64}(undef, nband, length(qx))
	for qxi in eachindex(qx)
		aimqdirects = [Vec3([qx[qxi] // Nqx, ky // Nqy, 0]) for ky in start[2]-1:stop[2]]
		aimq_qindex = Vector{Int}(undef, Nqy + 1)
		for (qi, qv) in enumerate(aimqdirects)
			aimq_qindex[qi] = findfirst(q -> all(isinteger, q - qv), redqgrid.kdirect)
		end
		G = Int.(aimqdirects[1] - redqgrid.kdirect[aimq_qindex[1]])
		nnkpts[5, 1] = aimq_qindex[1]
		nnkpts[6:8, 1] = G
		for i in 2:Nqy
			G = Int.(aimqdirects[i] - redqgrid.kdirect[aimq_qindex[i]])
			nnkpts[1, i-1] = aimq_qindex[i]
			nnkpts[2:4, i-1] = G
			nnkpts[5, i] = aimq_qindex[i]
			nnkpts[6:8, i] = G
		end
		G = Int.(aimqdirects[end] - redqgrid.kdirect[aimq_qindex[end]])
		nnkpts[1, end] = aimq_qindex[end]
		nnkpts[2:4, end] = G

		M = _Wilsonloop_mmn_BSE_NP(qgrid, bse, BSEband, bandindex, nnkpts)

		W = 1
		for i in 1:Nqy
			F = svd(M[:, :, i])
			W = (F.U * F.Vt) * W
		end

		θ_qx[:, qxi] = eigvals!(Hermitian(W))
	end
	wcc_qx = θ_qx ./ 2π

	return wcc_qx, wcc_qy
end


function _Wilsonloop_mmn_BSE_NP(qgrid, bse, BSEband, bandindex, nnkpts)

	qgrid = reducible_kgrid(qgrid)
	kgrid = reducible_kgrid(bse.kgrid)

	uband = BAND(kgrid, bse.TB, bse.TB.orb_location; vector = true)

	#A[vck] = A[i]
	nv = length(bse.v)
	nc = length(bse.c)
	Nk = length(kgrid)
	vck2i = Array{typeof(length(bse.vckindex))}(undef, nv, nc, Nk)
	for (i, (v, c, k)) in enumerate(bse.vckindex)
		vck2i[v, c, k] = i
	end

	vvcck = CartesianIndices((nv, nv, nc, nc, Nk))


	nband = length(bandindex)
	norb_ele = length(bse.TB.orb_location)
	nqb = size(nnkpts, 2)
	tasks = Array{Task}(undef, nband, nband, nqb)

	for i in 1:nqb
		q_qindex = nnkpts[1, i]
		q = qgrid.kdirect[q_qindex] + nnkpts[2:4, i]
		q_plus_b_qindex = nnkpts[5, i]
		qb = qgrid.kdirect[q_plus_b_qindex] + nnkpts[6:8, i]
		b2 = (qb - q) .// 2

		k_minus_b2_kindex = Vector{Int}(undef, Nk)
		k_plus_q_kindex = Vector{Int}(undef, Nk)
		kb2_plus_qb_kindex = Vector{Int}(undef, Nk)
		ΔG = Vector{Vec3{Int}}(undef, Nk)
		for (ki, kv) in enumerate(kgrid.kdirect)
			k_minus_b2 = kv - b2
			k_minus_b2_kindex[ki] = findfirst(k -> all(isinteger, k - k_minus_b2), kgrid.kdirect)

			k_plus_q = kv + q
			k_plus_q_kindex[ki] = findfirst(k -> all(isinteger, k - k_plus_q), kgrid.kdirect)
			k_plus_q_G = k_plus_q - kgrid.kdirect[k_plus_q_kindex[ki]]

			kb2_plus_qb = kgrid.kdirect[k_minus_b2_kindex[ki]] + qb
			kb2_plus_qb_kindex[ki] = findfirst(k -> all(isinteger, k - kb2_plus_qb), kgrid.kdirect)
			kb2_plus_qb_G = kb2_plus_qb - kgrid.kdirect[kb2_plus_qb_kindex[ki]]
			ΔG[ki] = k_plus_q_G - kb2_plus_qb_G
		end

		for α in 1:nband, β in 1:nband
			tasks[α, β, i] = Threads.@spawn sum(vvcck) do I
				(v, v′, c, c′, k) = Tuple(I)
				conj(BSEband[q_qindex].vectors[vck2i[v, c, k], bandindex[α]]) *
				BSEband[q_plus_b_qindex].vectors[vck2i[v′, c′, k_minus_b2_kindex[k]], bandindex[β]] *
				(sum(ii -> cis(2π * (ΔG[k] ⋅ bse.TB.orb_location[ii])) * conj(uband[k_plus_q_kindex[k]].vectors[ii, bse.c[c]]) * uband[kb2_plus_qb_kindex[k]].vectors[ii, bse.c[c′]], Base.OneTo(norb_ele))) *
				(uband[k_minus_b2_kindex[k]].vectors[:, bse.vindex[v′]] ⋅ uband[k].vectors[:, bse.vindex[v]])
			end
		end
	end

	M = fetch.(tasks)

	return M
end
