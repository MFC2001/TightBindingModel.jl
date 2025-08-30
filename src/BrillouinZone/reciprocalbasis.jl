function reciprocalbasis(lattvec::AbstractMatrix{<:Real}, periodicity)

	p = count(periodicity .== "p")
	if p == 3
		a₁ = lattvec[:, 1]
		a₂ = lattvec[:, 2]
		a₃ = lattvec[:, 3]

		V = a₁ ⋅ (a₂ × a₃)
		b₁ = (a₂ × a₃) * 2π / V
		b₂ = (a₃ × a₁) * 2π / V
		b₃ = (a₁ × a₂) * 2π / V
		U = [b₁ b₂ b₃]

	elseif p == 2
		if periodicity[1] == "np"
			a₂ = lattvec[:, 2]
			a₃ = lattvec[:, 3]
			a₁ = a₂ × a₃

			V = a₁ ⋅ (a₂ × a₃)
			b₂ = (a₃ × a₁) * 2π / V
			b₃ = (a₁ × a₂) * 2π / V
			U = [zeros(3) b₂ b₃]

		elseif periodicity[2] == "np"
			a₁ = lattvec[:, 1]
			a₃ = lattvec[:, 3]
			#Watch out!
			a₂ = a₃ × a₁

			V = a₁ ⋅ (a₂ × a₃)
			b₁ = (a₂ × a₃) * 2π / V
			b₃ = (a₁ × a₂) * 2π / V
			U = [b₁ zeros(3) b₃]

		elseif periodicity[3] == "np"
			a₁ = lattvec[:, 1]
			a₂ = lattvec[:, 2]
			a₃ = a₁ × a₂

			V = a₁ ⋅ (a₂ × a₃)
			b₁ = (a₂ × a₃) * 2π / V
			b₂ = (a₃ × a₁) * 2π / V
			U = [b₁ b₂ zeros(3)]

		else
			error("Wrong periodicity from kline.")
		end

	elseif p == 1
		if periodicity[1] == "p"
			a = lattvec[:, 1]
			b = 2π / sum(abs2, a) * a
			U = [b zeros(3, 2)]

		elseif periodicity[2] == "p"
			a = lattvec[:, 2]
			b = 2π / sum(abs2, a) * a
			U = [zeros(3) b zeros(3)]

		elseif periodicity[3] == "p"
			a = lattvec[:, 3]
			b = 2π / sum(abs2, a) * a
			U = [zeros(3, 2) b]

		else
			error("Wrong periodicity from kline.")
		end

	elseif p == 0
		U = zeros(Int, 3, 3)
	else
		error("Wrong periodicity from kline.")
	end

	return U
end

