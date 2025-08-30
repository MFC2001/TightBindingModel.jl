export Gaussian1D, Gaussian3D

function Gaussian1D(; r₀ = 0, σ = 1)
	A = 1 / (√(2π) * σ)

	GF = function (x)
		A * exp(-(x -  r₀)^2 / (2 * σ^2))
	end
	return GF
end

function Gaussian3D(; r₀ = [0, 0, 0], σ = [1, 1, 1])
	A = 1 / ((2π)^(3 / 2) * prod(σ))

	GF = function (r)
		return A * exp(-sum(((r - r₀) ./ σ) .^ 2) / 2)
	end

	return GF
end
