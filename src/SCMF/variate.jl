export Δ2x, x2Δ
function Δ2x(Δ, norb, Nk)
	x = Vector{real(eltype(Δ[1]))}(undef, Nk * norb * norb)

	t = 1
	for k in 1:Nk
		for j in 2:norb, i in 1:j-1
			(x[t], x[t+1]) = reim(Δ[k][i, j])
			t += 2
		end
		for i in 1:norb
			x[t] = real(Δ[k][i, i])
			t += 1
		end
	end

	return x
end
function x2Δ(x, norb, Nk)
	type = complex(eltype(x))
	Δ = Vector{Matrix{type}}(undef, Nk)

	t = 1
	for k in 1:Nk
		Δk = Matrix{type}(undef, norb, norb)
		for j in 2:norb, i in 1:j-1
			Δk[i, j] = complex(x[t], x[t+1])
			Δk[j, i] = complex(x[t], -x[t+1])
			t += 2
		end
		for i in 1:norb
			Δk[i, i] = complex(x[t])
			t += 1
		end
		Δ[k] = Δk
	end

	return Δ
end
