
function _mmn_BSE(qgrid, bse, BSEband, bandindex, nnkpts)
	if _mmn_BSE_is_double_kgrid(qgrid, bse.kgrid, bse.TB.period)
		return _mmn_BSE_double_kgrid(qgrid, bse, BSEband, bandindex, nnkpts)
	else
		return _mmn_BSE_random_kgrid(qgrid, bse, BSEband, bandindex, nnkpts)
	end
end

function _mmn_BSE_double_kgrid(qgrid, bse, BSEband, bandindex, nnkpts)

	#Find periodic direction
	Tp = _mmn_BSE_double_kgrid_findp(bse.TB.period)

	kgrid = bse.kgrid
	vckmap = bse.vckmap
	(nv, nc, nk) = size(vckmap)

	nq = length(BSEband)
	nband = length(bandindex)
	norb_ele = length(bse.TB.orb_location)
	nqb = size(nnkpts, 2)


	k_plus_q_kindex = Matrix{Int}(undef, nk, nqb)
	k_plus_q_plus_b2_kindex = Matrix{Int}(undef, nk, nqb)
	k_minus_b2_kindex = Matrix{Int}(undef, nk, nqb)
	b2 = Vector{Vec3{Rational{Int}}}(undef, nqb)
	Threads.@threads for i in 1:nqb
		q_qindex = nnkpts[1, i]
		q = qgrid.kdirect[q_qindex]
		q_plus_b_qindex = nnkpts[2, i]
		q_plus_b = qgrid.kdirect[q_plus_b_qindex] + nnkpts[3:5, i]
		b2[i] = (q_plus_b - q) .// 2
		b2_p = [0.0, 0.0, 0.0]
		b2_p[Tp] .= b2[Tp]
		for (ki, k) in enumerate(kgrid)
			# k = kgrid[ki]

			k_plus_q = k + q
			k_plus_q_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_plus_q), kgrid)

			k_plus_q_plus_b2 = k + q_plus_b - b2_p
			k_plus_q_plus_b2_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_plus_q_plus_b2), kgrid)

			k_minus_b2 = k - b2_p
			k_minus_b2_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_minus_b2), kgrid)
		end
	end

	#Calculate Mc and Mv
	eb2orb = Matrix{ComplexF64}(undef, norb_ele, nqb)
	Threads.@threads for i in 1:nqb
		eb2orb[:, i] = map(x -> cis(-2π * (b2[i] ⋅ x)), bse.TB.orb_location)
	end
	band = bse.bandk
	#Mc
	tasks = Array{Task}(undef, nc, nc, nk, nqb)
	for i in 1:nqb, k in 1:nk, c′ in 1:nc, c in 1:nc
		tasks[c, c′, k, i] = Threads.@spawn sum(Base.OneTo(norb_ele)) do iorb
			eb2orb[iorb, i] * conj(band[k_plus_q_kindex[k, i]].vectors[iorb, vckmap.idx2c[c]]) * band[k_plus_q_plus_b2_kindex[k, i]].vectors[iorb, vckmap.idx2c[c′]]
		end
	end
	Mc = fetch.(tasks)
	#Mv
	tasks = Array{Task}(undef, nv, nv, nk, nqb)
	for i in 1:nqb, k in 1:nk, v′ in 1:nv, v in 1:nv
		tasks[v, v′, k, i] = Threads.@spawn sum(Base.OneTo(norb_ele)) do iorb
			eb2orb[iorb, i] * conj(band[k_minus_b2_kindex[k, i]].vectors[iorb, vckmap.idx2v[v′]]) * band[k].vectors[iorb, vckmap.idx2v[v]]
		end
	end
	Mv = fetch.(tasks)

	#Calculate M
	nq = length(BSEband)
	Avck = Array{ComplexF64}(undef, nv, nc, nk, nband, nq)
	for q in 1:nq, α in 1:nband, k in 1:nk, c in 1:nc, v in 1:nv
		Avck[v, c, k, α, q] = BSEband[q].vectors[vckmap.vck2idx[v, c, k], bandindex[α]]
	end

	vvcck = Tuple.(CartesianIndices((nv, nv, nc, nc, nk)))

	tasks = Array{Task}(undef, nband, nband, nqb)
	for i in 1:nqb, β in 1:nband, α in 1:nband
		tasks[α, β, i] = Threads.@spawn sum(vvcck) do (v, v′, c, c′, k)
			conj(Avck[v, c, k, α, nnkpts[1, i]]) * Avck[v′, c′, k_minus_b2_kindex[k, i], β, nnkpts[2, i]] *
			Mc[c, c′, k, i] * Mv[v, v′, k, i]
		end
	end
	M = fetch.(tasks)

	return M
end
function _mmn_BSE_random_kgrid(qgrid, bse, BSEband, bandindex, nnkpts)

	kgrid = bse.kgrid
	vckmap = bse.vckmap
	(nv, nc, nk) = size(vckmap)

	nq = length(BSEband)
	nband = length(bandindex)
	norb_ele = length(bse.TB.orb_location)
	nqb = size(nnkpts, 2)


	k_plus_q_kindex = Matrix{Int}(undef, nk, nqb)
	k_plus_q_plus_b_kindex = Matrix{Int}(undef, nk, nqb)
	b2 = Vector{Vec3{Rational{Int}}}(undef, nqb)
	Threads.@threads for i in 1:nqb
		q_qindex = nnkpts[1, i]
		q = qgrid.kdirect[q_qindex]
		q_plus_b_qindex = nnkpts[2, i]
		q_plus_b = qgrid.kdirect[q_plus_b_qindex] + nnkpts[3:5, i]
		b2[i] = (q_plus_b - q) .// 2
		for (ki, k) in enumerate(kgrid)
			# k = kgrid[ki]

			k_plus_q = k + q
			k_plus_q_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_plus_q), kgrid)

			k_plus_q_plus_b = k + q_plus_b
			k_plus_q_plus_b_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_plus_q_plus_b), kgrid)
		end
	end

	#Calculate Pk
	R = gridindex(kgrid.kgrid_size)
	tasks = Array{Task}(undef, nk, nk, nqb)
	for i in 1:nqb, k′ in 1:nk, k in 1:nk
		tasks[k, k′, i] = Threads.@spawn begin
			kkb2 = kgrid[k′] - kgrid[k] + b2[i]
			sum(R -> cis(2π * (kkb2 ⋅ R)), R) / nk
		end
	end
	Pk = fetch.(tasks)

	#Calculate Mc and Mv
	eb2orb = Matrix{ComplexF64}(undef, norb_ele, nqb)
	Threads.@threads for i in 1:nqb
		eb2orb[:, i] = map(x -> cis(-2π * (b2[i] ⋅ x)), bse.TB.orb_location)
	end
	band = bse.bandk
	#Mc
	tasks = Array{Task}(undef, nc, nc, nk, nk, nqb)
	for i in 1:nqb, k′ in 1:nk, k in 1:nk, c′ in 1:nc, c in 1:nc
		tasks[c, c′, k, k′, i] = Threads.@spawn sum(Base.OneTo(norb_ele)) do iorb
			eb2orb[iorb, i] * conj(band[k_plus_q_kindex[k, i]].vectors[iorb, vckmap.idx2c[c]]) * band[k_plus_q_plus_b_kindex[k′, i]].vectors[iorb, vckmap.idx2c[c′]]
		end
	end
	Mc = fetch.(tasks)
	#Mv
	tasks = Array{Task}(undef, nv, nv, nk, nk, nqb)
	for i in 1:nqb, k′ in 1:nk, k in 1:nk, v′ in 1:nv, v in 1:nv
		tasks[v, v′, k, k′, i] = Threads.@spawn sum(Base.OneTo(norb_ele)) do iorb
			eb2orb[iorb, i] * conj(band[k′].vectors[iorb, vckmap.idx2v[v′]]) * band[k].vectors[iorb, vckmap.idx2v[v]]
		end
	end
	Mv = fetch.(tasks)

	nq = length(BSEband)
	Avck = Array{ComplexF64}(undef, nv, nc, nk, nband, nq)
	for q in 1:nq, α in 1:nband, k in 1:nk, c in 1:nc, v in 1:nv
		Avck[v, c, k, α, q] = BSEband[q].vectors[vckmap.vck2idx[v, c, k], bandindex[α]]
	end

	vvcckk = Tuple.(CartesianIndices((nv, nv, nc, nc, nk, nk)))

	tasks = Array{Task}(undef, nband, nband, nqb)
	for i in 1:nqb, β in 1:nband, α in 1:nband
		tasks[α, β, i] = Threads.@spawn sum(vvcckk) do (v, v′, c, c′, k, k′)
			conj(Avck[v, c, k, α, nnkpts[1, i]]) * Avck[v′, c′, k′, β, nnkpts[2, i]] *
			Pk[k, k′, i] * Mc[c, c′, k, k′, i] * Mv[v, v′, k, k′, i]
		end
	end
	M = fetch.(tasks)

	return M
end

function _mmn_BSE_is_double_kgrid(qgrid, kgrid, period)
	p = count(period)
	if p == 3
		return qgrid.kgrid_size == kgrid.kgrid_size .* 2
	elseif p == 2
		if !period[1]
			T = [2, 3]
		elseif !period[2]
			T = [1, 3]
		elseif !period[3]
			T = [1, 2]
		end
		return qgrid.kgrid_size[T] == kgrid.kgrid_size[T] .* 2
	elseif p == 1
		if period[1]
			return qgrid.kgrid_size[1] == kgrid.kgrid_size[1] .* 2
		elseif period[2]
			return qgrid.kgrid_size[2] == kgrid.kgrid_size[2] .* 2
		elseif period[3]
			return qgrid.kgrid_size[3] == kgrid.kgrid_size[3] .* 2
		end
	elseif p == 0
		return true
	end
end
function _mmn_BSE_double_kgrid_findp(period)
	p = count(period)
	if p == 3
		Tp = [1, 2, 3]
	elseif p == 2
		if !period[1]
			Tp = [2, 3]
		elseif !period[2]
			Tp = [1, 3]
		elseif !period[3]
			Tp = [1, 2]
		end
	elseif p == 1
		if period[1]
			Tp = [1]
		elseif period[2]
			Tp = [2]
		elseif period[3]
			Tp = [3]
		end
	elseif p == 0
		Tp = Int[]
	end
	return Tp
end
