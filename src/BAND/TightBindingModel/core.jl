
function BAND(kpoints, rh::HermitianReciprocalHoppings, ::Val{true})
	band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[i] = eigen!(rh(kpoints[i]))
	end
	return band
end
function BAND(kpoints, rh::HermitianReciprocalHoppings, ::Val{false})
	band = Matrix{Float64}(undef, numorb(rh), length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[:, i] = eigvals!(rh(kpoints[i]))
	end
	return band
end
function BAND(kpoints, rh::HermitianReciprocalHoppings, orblocat, ::Val{true})
	band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[i] = eigen!(rh(kpoints[i]), orblocat)
	end
	return band
end
function BAND(kpoints, rh::HermitianReciprocalHoppings, orblocat, ::Val{false})
	band = Matrix{Float64}(undef, numorb(rh), length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[:, i] = eigvals!(rh(kpoints[i]), orblocat)
	end
	return band
end
