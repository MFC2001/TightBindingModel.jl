function istaskrunnable(task::Task)
	return istaskstarted(task) && !istaskdone(task)
end



# 关于求解特征值的加速方案，一个可能的实现方式
# 但是，使用julia原生的包时，底层都有BLAS并行加速，至少目前该方案并没有达到明显缩短运行时间的目的。
function BSE_NP(
	bse::BSE,
	qpoints::AbstractVector{<:ReducedCoordinates{<:Real}};
	vector = false,
)

	kgrid = bse.kgrid
	nk = length(kgrid)
	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	addmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, +)
	minusmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, -)

	norb = numorb(bse.TB)
	U_screened_k = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Threads.@threads for k in 1:nk
		bse.U_screened!(U_screened_k[k], kgrid_Γ[k])
		U_screened_k[k] ./= nk
	end
	J_unscreened_k = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Threads.@threads for k in 1:nk
		bse.J_unscreened!(J_unscreened_k[k], kgrid_Γ[k])
		J_unscreened_k[k] ./= nk
	end

	U_unscreened_q = Matrix{ComplexF64}(undef, norb, norb)
	J_screened_kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	J_unscreened_kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]


	Nq = length(qpoints)
	tasks = Vector{Task}(undef, 2 * Nq)
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	for (qi, q) in enumerate(qpoints)

		# We can't calculate Γ directly.
		if norm(q) == 0
			p = bse.TB.period .== "p"
			if count(p) == 0
				q = [0.0, 0.0, 0.0]
			else
				I = findfirst(p)
				q = [0.0, 0.0, 0.0]
				q[I] = 1e-6
			end
		end

		bandkq = BAND(map(k -> k + q, kgrid), bse.TB; vector = true)

		task = Threads.@spawn begin
			bse.U_unscreened!(U_unscreened_q, q)
			U_unscreened_q ./= nk
		end
		Threads.@threads for k in 1:nk
			bse.J_screened!(J_screened_kq[k], kgrid_Γ[k] + q)
			J_screened_kq[k] ./= nk
		end
		Threads.@threads for k in 1:nk
			bse.J_unscreened!(J_unscreened_kq[k], kgrid_Γ[k] + q)
			J_unscreened_kq[k] ./= nk
		end
		wait(task)

		Kᵈ = CreatKᵈ_Q(bse.bandk, bandkq, U_screened_k, J_screened_kq[minusmap[1, 1]], J_screened_kq, addmap, minusmap)
		Kˣ = CreatKˣ_Q(bse.bandk, bandkq, U_unscreened_q, J_unscreened_k, J_unscreened_kq, addmap, minusmap)

		H_t = nothing
		H_s = nothing
		while true
			try
				H_t, H_s = BSE_NP_Hamilton(bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)
				break
			catch err
				if err isa OutOfMemoryError
					ntaskrunnable = count(istaskrunnable, tasks[1:2*(qi-1)])
					if ntaskrunnable > 0
						@info "qi = $qi \n ntaskrunnable = $ntaskrunnable \n Wait an eigen task done for enough memory!"
						while true
							if count(istaskrunnable, tasks[1:2*(qi-1)]) == ntaskrunnable
								sleep(5)
							else
								break
							end
						end
					else
						rethrow(err)
					end
				else
					rethrow(err)
				end
			end
		end

		tasks[2*qi-1] = Threads.@spawn begin
			BSEband_t[qi] = _eigsolve_Hmat(H_t)
		end
		tasks[2*qi] = Threads.@spawn begin
			BSEband_s[qi] = _eigsolve_Hmat(H_s)
		end
	end
	wait.(tasks)

	if !vector
		BSEband_t_vals
	end

	return BSEband_t, BSEband_s
end