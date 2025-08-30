export DOS

function DOS(energy::AbstractArray{<:Real}; η = 0.01, E = minimum(energy)-8*η:η/5:maximum(energy)+8*η)

	density = zeros(Float64, length(E))

	for ε in energy
		F = Gaussian1D(; b = ε, c = η)
		density .+= F.(E)
	end

	return E, density
end
