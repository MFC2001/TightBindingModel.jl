
struct RealSpaceWannierExciton <: RealSpaceExciton
	nq::Int
	nk::Int
	qgrid::RedKgrid
	kgrid::RedKgrid
	kq_kindex::Matrix{Int}
	WU::Array{ComplexF64} #∑_{vc} WUU* 
end
function (w::RealSpaceWannierExciton)(ei, Rₑ, hi, Rₕ, R = [0, 0, 0])
	w_iR = Vector{ComplexF64}(undef, w.nq)
	Threads.@threads for q in Base.OneTo(w.nq)
		w_iR[q] = sum(Base.OneTo(w.nk)) do k
			w.WU[k, q, ei, hi] * cis(2π * (w.kgrid[w.kq_kindex[k, q]] ⋅ Rₑ - w.kgrid[k] ⋅ Rₕ))
		end
		w_iR[q] *= cis(-2π * (w.qgrid[q] ⋅ R))
	end
	return sum(w_iR) / w.nk
end

function RealSpaceWannierExciton(qgrid::RedKgrid, bse::AbstractBSE, BSEband::AbstractVector{<:Eigen}, u_matrix::Array{<:Number, 3};
	u_matrix_opt::Array{<:Number, 3} = Array{Int, 3}(0, 0, 0), exclude_bands::AbstractVector{<:Integer} = Int[])

	BSEband_opt = disentangled_band(BSEband, u_matrix_opt, exclude_bands)
	wannier = wannierize(BSEband_opt, u_matrix)
	(_, nq, nw) = size(wannier)

	vckmap = bse.vckmap
	kgrid = bse.kgrid
	bandk = bse.bandk
	(nv, nc, nk) = size(vckmap)


	k_plus_q_kindex = Matrix{Int}(undef, nk, nq)
	Threads.@threads for qi in 1:nq
		q = qgrid[qi]
		for (ki, k) in enumerate(kgrid.kdirect)
			# k = kgrid.kdirect[ki]
			k_plus_q = k + q
			k_plus_q_kindex[ki, qi] = findfirst(k -> all(isinteger, k - k_plus_q), kgrid.kdirect)
		end
	end

	wannier_vck = Array{eltype(wannier)}(undef, nv, nc, nk, nq, nw)
	Threads.@threads for qw in CartesianIndices((nq, nw))
		for k in 1:nk, c in 1:nc, v in 1:nv
			wannier_vck[v, c, k, qw] = wannier[vckmap.vck2idx[v, c, k], qw]
		end
	end

	vc = [(v, c) for v in 1:nv, c in 1:nc]

	norb_ele = numorb(bse.TB)
	tasks = Array{Task}(undef, nk, nq, norb_ele, norb_ele)
	wannier_iR = Vector{RealSpaceWannierExciton}(undef, nw)
	for w in 1:nw
		for hi in 1:norb_ele, ei in 1:norb_ele, q in 1:nq, k in 1:nk
			tasks[k, q, ei, hi] = Threads.@spawn sum(vc) do (v, c)
				wannier_vck[v, c, k, q, w] *
				bandk[k_plus_q_kindex[k, q]].vectors[ei, vckmap.idx2c[c]] *
				conj(bandk[k].vectors[hi, vckmap.idx2v[v]])
			end
		end
		WU = fetch.(tasks)
		wannier_iR[w] = RealSpaceWannierExciton(nq, nk, qgrid, kgrid, k_plus_q_kindex, WU)
	end

	return wannier_iR
end



