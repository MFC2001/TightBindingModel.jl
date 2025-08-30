
struct RealSpaceBlochExciton <: RealSpaceExciton
	q::Vector{Float64}
	nk::Int
	inv_nk::Float64
	kgrid::RedKgrid
	AU::Array{ComplexF64} #∑_{vc} AUU* 
end
function (ψ::RealSpaceBlochExciton)(ei, Rₑ, hi, Rₕ)
	return ψ.inv_nk * sum(Base.OneTo(ψ.nk)) do k
		ψ.AU[k, ei, hi] * cis(2π * (ψ.kgrid[k] ⋅ (Rₑ - Rₕ) + q ⋅ Rₑ))
	end
end
function RealSpaceBlochExciton(bse::AbstractBSE, q, A)

	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)
	end

	vckmap = bse.vckmap
	kgrid = bse.kgrid

	bandkq = BAND(map(k -> k + q, kgrid), bse.TB; vector = true)
	_sum_wave_is_real!.(bandkq)


	(nv, nc, nk) = size(vckmap)

	Avck = Array{eltype(A)}(undef, nv, nc, nk)
	for k in 1:nk, c in 1:nc, v in 1:nv
		Avck[v, c, k] = A[vckmap.vck2idx[v, c, k]]
	end

	vc = [(v, c) for v in 1:nv, c in 1:nc]

	nk = length(kgrid)
	norb_ele = numorb(bse.TB)
	tasks = Array{Task}(undef, nk, norb_ele, norb_ele)
	for hi in 1:norb_ele, ei in 1:norb_ele, k in 1:nk
		tasks[k, ei, hi] = Threads.@spawn sum(vc) do (v, c)
			Avck[v, c, k] * bandkq[k].vectors[ei, vckmap.idx2c[c]] * conj(bse.bandk[k].vectors[hi, vckmap.idx2v[v]])
		end
	end
	AU = fetch.(tasks)

	return RealSpaceBlochExciton(q, nk, 1 / nk, kgrid, AU)
end
