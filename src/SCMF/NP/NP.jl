
export SCMF_NP
include("./order.jl")

"""

There are still some questions in this model.
"""
function SCMF_NP(
	symTB::SymTightBindModel,
	kgrid::Union{RedKgrid, IrredKgrid},
	Distribution::Function,
	Uhr::HR,
	Jhr::Union{HR, Nothing};
	Δ₀ = nothing,
	Δeps = 1e-10,
	maxittr = 100,
)


	nΔ = length(symTB.value)
	R = gridindex(kgrid.kgrid_size)


	Δindex = Vector{Tuple{Int, Int, Int}}(undef, nΔ)
	for (i, hop) in enumerate(symTB.irredhop)
		Ri = findfirst(x -> all(iszero, x - hop.R), R)
		Δindex[i] = (hop.i, hop.j, Ri)
	end

	ViR = CreatViR_NP(R, Uhr, Jhr)

	K = SCMF_NP_kernal(ViR, Δindex, R, numorb(symTB))

	DM = DensityMatrix_R(kgrid, Distribution)

	h₀ = deepcopy(symTB.value)
	if isnothing(Δ₀)
		V = rand(Float64, nΔ) .% 1
		Δ₀ = h₀ .* (1 .+ V)
		# Δ₀ = symhr.value
	else
		Δ₀ = convert.(eltype(h₀), Δ₀)
	end

	ΓR = findfirst(iszero, R)


	T = deepcopy(h₀)
	SCMFsymTB = deepcopy(symTB)
	F! = function (out, x)

		SCMFsymTB.value .= h₀ + x
		refresh!(SCMFsymTB)

		band = BAND(kgrid, SCMFsymTB; vector = true)
		ρ = DM(band)
		Threads.@threads for i in 1:nΔ
			T[i] = real(sum(K[i, :, :, :] .* ρ))
		end

		# T[1] = 0

		out .= T - x

		return nothing
	end

	result = NLsolve.nlsolve(F!, Δ₀, iterations = maxittr, method = :newton, ftol = Δeps, show_trace = true)

	SCMFsymTB.value .= h₀ + result.zero
	refresh!(SCMFsymTB)

	return SCMFsymTB
end

