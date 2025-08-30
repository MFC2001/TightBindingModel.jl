export CreatViR_NP
function CreatViR_NP(R::AbstractVector, Uhr::HR{TU}, Jhr::Union{HR, Nothing})::Function where {TU}

	norb = numorb(Uhr)

	minusmap = gridmap(R, -)

	UR = interaction_hr2array(Uhr, R)

	onsite = Vector{TU}(undef, norb)
	ΓRindex = minusmap[1, 1]
	for i in eachindex(onsite)
		onsite[i] = UR[i, i, ΓRindex]
		UR[i, i, ΓRindex] = 0
	end


	if isnothing(Jhr)
		ViR = function (i, R₁, j, R₂, k, R₃, l, R₄; spin = false)
			return (spin ? 0 : kronecker_product((i, l), (i, j), (i, k), (R₁, R₄), (R₁, R₂), (R₁, R₃)) * onsite[i]) +
				   kronecker_product((i, l), (R₁, R₄), (j, k), (R₂, R₃)) * UR[i, j, minusmap[R₂, R₁]]
		end
	else
		JR = interaction_hr2array(Jhr, R)
		ViR = function (i, R₁, j, R₂, k, R₃, l, R₄; spin = false)
			return (spin ? 0 : kronecker_product((i, l), (i, j), (i, k), (R₁, R₄), (R₁, R₂), (R₁, R₃)) * onsite[i]) +
				   kronecker_product((i, l), (R₁, R₄), (j, k), (R₂, R₃)) * UR[i, j, minusmap[R₂, R₁]] +
				   kronecker_product((i, k), (R₁, R₃), (j, l), (R₂, R₄)) * JR[i, j, minusmap[R₂, R₁]] +
				   kronecker_product((i, j), (R₁, R₂), (k, l), (R₃, R₄)) * JR[i, j, minusmap[R₂, R₁]]
		end
	end

	return ViR
end

function interaction_hr2array(hr::HR{T}, R) where {T}
	norb = numorb(hr)
	array = zeros(T, norb, norb, length(R))
	for (i, v) in enumerate(hr.value)
		Ri = findfirst(x -> x[1] == hr.path[i, 1] && x[2] == hr.path[i, 2] && x[3] == hr.path[i, 3], R)
		array[hr.path[i, 4], hr.path[i, 5], Ri] = v
	end
	return array
end
