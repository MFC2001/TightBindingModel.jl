export BAND

"""
	BAND(kpoints, hr or rh; vector::Bool = false)
	BAND(kpoints, hr or rh, orblocat; vector::Bool = false)
"""


function BAND(kgrid::MonkhorstPack, args...; vector::Bool = false)
	return BAND(RedKgrid(kgrid), args...; vector)
end
function BAND(kpoints::AbstractBrillouinZone, args...; vector::Bool = false)
	return BAND(kpoints.kdirect, args...; vector)
end

function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, TB::AbstractTightBindModel; vector::Bool = false)
	return BAND(kpoints, TB.H, Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, TB::AbstractTightBindModel, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)
	return BAND(kpoints, TB.H, orblocat, Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, hr::HR; vector::Bool = false)
	return BAND(kpoints, HermitianReciprocalHoppings(hr), Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, hr::HR, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)
	return BAND(kpoints, HermitianReciprocalHoppings(hr), orblocat, Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, rh::HermitianReciprocalHoppings; vector::Bool = false)
	return BAND(kpoints, rh, Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, rh::HermitianReciprocalHoppings, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)
	return BAND(kpoints, rh, orblocat, Val(vector))
end
function BAND(kpoints, rh, ::Val{true})
	band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[i] = eigen!(rh(kpoints[i]))
	end
	return band
end
function BAND(kpoints, rh, ::Val{false})
	band = Matrix{Float64}(undef, numorb(hr), length(kdirects))
	Threads.@threads for i in eachindex(kdirects)
		band[:, i] = eigvals!(rh(kpoints[i]))
	end
	return band
end
function BAND(kpoints, rh, orblocat, ::Val{true})
	band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[i] = eigen!(rh(kpoints[i]), orblocat)
	end
	return band
end
function BAND(kpoints, rh, orblocat, ::Val{false})
	band = Matrix{Float64}(undef, numorb(hr), length(kdirects))
	Threads.@threads for i in eachindex(kdirects)
		band[:, i] = eigvals!(rh(kpoints[i]), orblocat)
	end
	return band
end













# function BAND_real(args...; vector = false)

# 	if vector
# 		H = Hamilton([0, 0, 0], args...)
# 		band = [eigen!(H)]
# 	else
# 		band = Matrix{Float64}(undef, numorb(args[1]), 1)
# 		H = Hamilton([0, 0, 0], args...)
# 		band[:, 1] = eigvals!(H)
# 	end

# 	return band
# end





# function BAND(Hk::AbstractVector{<:AbstractMatrix{<:Complex}}; vector = false)
# 	Nk = length(Hk)

# 	if vector
# 		band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nk)
# 		Threads.@threads for (i, H) in enumerate(Hk)
# 			band[i] = eigen!(Hermitian(H))
# 		end
# 	else
# 		band = Matrix{Float64}(undef, size(Hk[1], 1), Nk)
# 		Threads.@threads for (i, H) in enumerate(Hk)
# 			band[:, i] = eigvals!(Hermitian(H))
# 		end
# 	end

# 	return band
# end
# function BAND(Hk::AbstractVector{<:AbstractMatrix{<:Real}}; vector = false)
# 	Nk = length(Hk)

# 	if vector
# 		band = Vector{Eigen{Float64, Float64, Matrix{Float64}, Vector{Float64}}}(undef, Nk)
# 		Threads.@threads for (i, H) in enumerate(Hk)
# 			band[i] = eigen!(Hermitian(H))
# 		end
# 	else
# 		band = Matrix{Float64}(undef, size(Hk[1], 1), Nk)
# 		Threads.@threads for (i, H) in enumerate(Hk)
# 			band[:, i] = eigvals!(Hermitian(H))
# 		end
# 	end

# 	return band
# end
# function BAND(Hk::AbstractArray{<:Complex, 3}; vector = false)

# 	Nk = size(Hk, 3)

# 	if vector
# 		band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nk)
# 		Threads.@threads for i in 1:Nk
# 			band[i] = eigen!(Hermitian(Hk[:, :, i]))
# 		end
# 	else
# 		band = Matrix{Float64}(undef, size(Hk, 1), Nk)
# 		Threads.@threads for i in 1:Nk
# 			band[:, i] = eigvals!(Hermitian(Hk[:, :, i]))
# 		end
# 	end

# 	return band
# end
# function BAND(Hk::AbstractArray{<:Real, 3}; vector = false)

# 	Nk = size(Hk, 3)

# 	if vector
# 		band = Vector{Eigen{Float64, Float64, Matrix{Float64}, Vector{Float64}}}(undef, Nk)
# 		Threads.@threads for i in 1:Nk
# 			band[i] = eigen!(Hermitian(Hk[:, :, i]))
# 		end
# 	else
# 		band = Matrix{Float64}(undef, size(Hk, 1), Nk)
# 		Threads.@threads for i in 1:Nk
# 			band[:, i] = eigvals!(Hermitian(Hk[:, :, i]))
# 		end
# 	end

# 	return band
# end
