


function SCMF_kernal(H₀, order, Δ₀, Δeps, maxiter)

	Nk = length(H₀)
	norb = size(H₀[1], 1)

	band = Vector{Eigen}(undef, Nk)

	F! = function (out, x)

		Δx = x2Δ(x, norb, Nk)
		Threads.@threads for k in eachindex(band)
			band[k] = eigen!(H₀[k] + Hermitian(Δx[k]))
		end
		Δ = order(band)
		out .= x - Δ2x(Δ, norb, Nk)

		return nothing
	end

	result = NLsolve.nlsolve(F!, Δ2x(Δ₀, norb, Nk), iterations = maxiter, method = :newton, ftol = Δeps, show_trace = true)

	Δ = x2Δ(result.zero, norb, Nk)

	return [H₀[k] + Δ[k] for k in eachindex(H₀)]
end
