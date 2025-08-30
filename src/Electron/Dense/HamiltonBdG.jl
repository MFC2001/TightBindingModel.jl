
function HamiltonSC(hr::HR, Δ::Number)

	if Δ == 0
		H = Hamilton(hr)
	else

		He = Hamilton(hr) / 2
		Hh = -transpose(He)

		#Require the order of basis.
		norb = hr.norb


		iσy = [0 Δ; -Δ 0] / 2
		IΔ = blockdiagm(0 => repeat([iσy], norb ÷ 2))
		IΔ = Array(IΔ)

		H = BlockArray(undef_blocks, Matrix{ComplexF64}, [norb, norb], [norb, norb])
		H[Block(1, 1)] = He
		H[Block(2, 2)] = Hh
		H[Block(1, 2)] = IΔ
		H[Block(2, 1)] = IΔ'

		H = Array(H)
	end

	return Hermitian(H, :U)
end