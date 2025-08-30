
export smooth_AngleSingular

function smooth_AngleSingular(theta::AbstractVector{<:Real}, min::Real = 1)

	n = length(theta)
	dt = diff(theta)

	I1 = findall(dt .> min)
	I2 = findall(dt .< -min)

	for i in I1
		theta[i+1:end] .-= 2π
	end

	for i in I2
		theta[i+1:end] .+= 2π
	end

	return theta
end





