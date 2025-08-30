function _is_same_wave(ψ₁::AbstractVector{<:Number}, ψ₂::AbstractVector{<:Number}; atol = 1e-4)
	C = 1 / length(ψ₁)
	I = findfirst(x -> abs(x) >= C, ψ₁)
	T = ψ₁[I] / ψ₂[I]
	return all(i -> isapprox(ψ₁[i], ψ₂[i] * T; atol), eachindex(ψ₁))
end
