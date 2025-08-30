
function _amn_BSE(BSEband, bandindex, W)

	(_, nq, nw) = size(W)
	nband = length(bandindex)

	tasks = Array{Task}(undef, nband, nw, nq)
	for I in CartesianIndices(tasks)
		(α, wi, qi) = Tuple(I)
		tasks[I] = Threads.@spawn BSEband[qi].vectors[:, bandindex[α]] ⋅ W[:, qi, wi]
	end
	A = fetch.(tasks)

	return A
end

function _guess_iR_to_vck(qgrid, bse, guess::AbstractVector)

	#make sure that [0,0,0] exists.
	R = gridindex(bse.kgrid.kgrid_size)

	norb = length(bse.TB.orb_location)
	nR = length(R)
	nw = length(guess)

	#TODO 这里的类型需要考虑，目前假定输入实函数。
	WR = Array{Float64}(undef, norb, nR, norb, nR, nw)
	Threads.@threads for I in CartesianIndices((norb, nR, nw))
		(hi, Rₕ, wi) = Tuple(I)
		for Rₑ in Base.OneTo(nR), ei in Base.OneTo(norb)
			WR[ei, Rₑ, hi, Rₕ, wi] = guess[wi](ei, R[Rₑ], hi, R[Rₕ])
		end
	end

	# WR = Array{Float64}(undef, norb, NR, norb, NR, nw)
	# tasks = Array{Task}(undef, norb, NR, norb, nw)
	# for (wi, wf) in enumerate(guess), ei in Base.OneTo(norb), R_e in Base.OneTo(NR), hi in Base.OneTo(norb)
	# 	tasks[ei, R_e, hi, wi] = Threads.@spawn for R_h in Base.OneTo(NR)
	# 		WR[ei, R_e, hi, R_h, wi] = wf(ei, R[R_e], hi, R[R_h])
	# 	end
	# end
	# wait.(tasks)


	#Get main WR
	WR_main = Vector{Vector{eltype(WR)}}(undef, nw)
	WR_main_idx = Vector{Vector{Tuple{Int, Int, Int, Int}}}(undef, nw)
	Threads.@threads for wi in Base.OneTo(nw)
		WR_abs2 = map(abs2, WR[:, :, :, :, wi])
		WR_max = maximum(WR_abs2)
		WR_max_100 = WR_max / 100
		I = findall(x -> x > WR_max_100, WR_abs2)
		WR_main[wi] = WR[I, wi]
		WR_main_idx[wi] = Tuple.(I)
	end


	kgrid = bse.kgrid

	nq = length(qgrid)
	nk = length(kgrid)

	k_plus_q_kindex = Matrix{Int}(undef, nk, nq)
	Threads.@threads for qi in eachindex(qgrid)
		for ki in 1:nk
			k_plus_q = kgrid.kdirect[ki] + qgrid[qi]
			k_plus_q_kindex[ki, qi] = findfirst(k -> all(isinteger, k - k_plus_q), kgrid.kdirect)
		end
	end


	nR_inv = 1 / nR
	band = bse.bandk
	tasks = Array{Task}(undef, length(bse.vckmap), nq, nw)
	for wi in 1:nw, q in eachindex(qgrid), (i, (v, c, k)) in enumerate(bse.vckmap.idx2vck)
		tasks[i, q, wi] = Threads.@spawn begin
			kq = k_plus_q_kindex[k, q]
			U = conj.(band[kq].vectors[:, c]) * transpose(band[k].vectors[:, v])
			k_frac = kgrid.kdirect[k]
			kq_frac = kgrid.kdirect[kq]
			nR_inv * sum(eachindex(WR_main[wi])) do ii
				(ei, Rₑ, hi, Rₕ) = WR_main_idx[wi][ii]
				return WR_main[wi][ii] * U[ei, hi] * cis(2π * (k_frac ⋅ R[Rₕ] - kq_frac ⋅ R[Rₑ]))
			end
		end
	end
	W = fetch.(tasks)

	#TODO 应该在投影至vck基后归一化，对WR有两次筛选，一次是vc基，一次是q少于k。
	Threads.@threads for wi in 1:nw
		N = √sum(abs2, W[:, :, wi])
		W[:, :, wi] ./= N
	end

	return W
end
