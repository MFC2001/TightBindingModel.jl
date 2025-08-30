
function Exciton_NP_wannier_guess_Gauss(TB, center::AbstractVector, σR::AbstractVector, σr::AbstractVector = σR)

	if eltype(center) <: ReducedCoordinates{<:Real}
		center = map(x -> TB.lattice * x, center)
	elseif eltype(center) <: CartesianCoordinates{<:Real}
	else
		println("Make sure the center is a vector of CartesianCoordinates.")
	end

	guess = Vector{Function}(undef, length(center))

	for wi in eachindex(center)
		gaussian_R = Gaussian3D(; r₀ = center[wi], σ = σR[wi])
		gaussian_r = Gaussian3D(; r₀ = [0, 0, 0], σ = σr[wi])
		guess[wi] = function (ei, R_e, hi, R_h)
			r_e = TB.lattice * (TB.orb_location[ei] + R_e)
			r_h = TB.lattice * (TB.orb_location[hi] + R_h)
			return gaussian_R((r_e + r_h) / 2) * gaussian_r(r_e - r_h)
		end
	end

	return guess
end

function exciton_SP_wannier_guess_Gauss(TB, center::AbstractVector, σR::AbstractVector, σr::AbstractVector = σR)

	if eltype(center) <: ReducedCoordinates{<:Real}
		center = map(x -> TB.lattice * x, center)
	elseif eltype(center) <: CartesianCoordinates{<:Real}
	else
		println("Make sure the center is a vector of CartesianCoordinates.")
	end

	guess = Vector{Function}(undef, length(center))

	for wi in eachindex(center)
		gaussian_R = Gaussian3D(; r₀ = center[wi], σ = σR[wi])
		gaussian_r = Gaussian3D(; r₀ = [0, 0, 0], σ = σr[wi])
		guess[wi] = function (ei, R_e, hi, R_h)
			r_e = TB.lattice * (TB.orb_location[ei] + R_e)
			r_h = TB.lattice * (TB.orb_location[hi] + R_h)
			return gaussian_R((r_e + r_h) / 2) * gaussian_r(r_e - r_h)
		end
	end

	return guess
end