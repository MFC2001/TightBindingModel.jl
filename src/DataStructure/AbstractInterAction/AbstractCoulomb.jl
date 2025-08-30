
abstract type AbstractCoulomb <: AbstractInterAction end
abstract type ReciprocalCoulomb <: AbstractCoulomb end
abstract type RealCoulomb <: AbstractCoulomb end

#1e-19
const qₑ = 1.602176634
#1e-12
const ϵ₀ = 8.854187817
#This term equal to e^2/4πϵ₀ with the unit of r is Å.
#Coulomb potential is CoulombScale/r, the unit of r is Å, potential Energy unit is eV.
const CoulombScale = qₑ * 1e3 / (4 * π * ϵ₀)


struct RealInverseR <: RealCoulomb
	ϵ::Float64
	CoulombScaleϵ::Float64
end
function (v::RealInverseR)(r::AbstractVector{<:Real})
	r = norm(r)
	return r == 0 ? 0 : v.CoulombScaleϵ / r
end
function (v::RealInverseR)(r::Real)
	return r == 0 ? 0 : v.CoulombScaleϵ / r
end
function RealInverseR(; ϵ::Real = 1)
	CoulombScaleϵ = CoulombScale / ϵ
	return RealInverseR(ϵ, CoulombScaleϵ)
end

struct RealGauss <: RealCoulomb
	ϵ::Float64
	α::Float64
	CoulombScaleϵ::Float64
	V₀::Float64
end
function (v::RealGauss)(r::AbstractVector{<:Real})
	r = norm(r)
	return r == 0 ? v.V₀ : v.CoulombScaleϵ * erf(v.α * r) / r
end
function (v::RealGauss)(r::Real)
	return r == 0 ? v.V₀ : v.CoulombScaleϵ * erf(v.α * r) / r
end
function RealGauss(; ϵ::Real = 1, α::Real = 1)
	CoulombScaleϵ = CoulombScale / ϵ
	V₀ = CoulombScaleϵ * 2 * α / √π
	return RealGauss(ϵ, α, CoulombScaleϵ, V₀)
end


struct ReciprocalGauss3D <: ReciprocalCoulomb
	ϵ::Float64
	α::Float64
	Ω::Float64
	CoulombScaleϵ::Float64
	α²4::Float64
end
function (v::ReciprocalGauss3D)(k::AbstractVector{<:Real})::Float64
	k² = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]
	return k² == 0 ? 0.0 : v.CoulombScaleϵ * exp(-k² / v.α²4) / k²
end
function (v::ReciprocalGauss3D)(k::Real)::Float64
	if k == 0
		return 0.0
	else
		k² = k * k
		return v.CoulombScaleϵ * exp(-k² / v.α²4) / k²
	end
end
function (v::ReciprocalGauss3D)(sym::Symbol, nk::Integer)::Float64
	if sym == :head
		NΩ = nk * v.Ω
		qsz = (6 * π^2 / NΩ)^(1 // 3)
		return v.CoulombScaleϵ * NΩ * v.α * erf(qsz / (2 * v.α)) / (2 * π^(3 // 2))
	else
		error("Wrong method for a ReciprocalGauss3D.")
	end
end
function ReciprocalGauss3D(; ϵ::Real = 1, α::Real = 1, Ω::Real = 1)
	CoulombScaleϵ = CoulombScale * 4π / (Ω * ϵ)
	α²4 = 4 * α^2
	return ReciprocalGauss3D(ϵ, α, Ω, CoulombScaleϵ, α²4)
end


struct ReciprocalGauss2D <: ReciprocalCoulomb
	ϵ::Float64
	α::Float64
	S::Float64
	CoulombScaleϵ::Float64
	α2::Float64
end
function (v::ReciprocalGauss2D)(k::AbstractVector{<:Real}, z::Real)::Float64
	k² = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]

	if k² == 0
		return 0.0
	else
		k = √k²
		z = abs(z)

		αz = v.α * z
		ekz = exp(k * z)
		k2α = k / (v.α2)

		return v.CoulombScaleϵ * (erfc(k2α + αz) * ekz + erfc(k2α - αz) / ekz) / k
	end
end
function (v::ReciprocalGauss2D)(k::Real, z::Real)::Float64
	if k == 0
		return 0.0
	else
		z = abs(z)

		αz = v.α * z
		ekz = exp(k * z)
		k2α = k / (v.α2)

		return v.CoulombScaleϵ * (erfc(k2α + αz) * ekz + erfc(k2α - αz) / ekz) / k
	end
end
function (v::ReciprocalGauss2D)(sym::Symbol, nk::Integer, z::Real)::Float64
	if sym == :head
		NS = nk * v.S
		qsz = √(4 * π / NS)

		z = abs(z)

		qsz2α = qsz / (v.α2)
		if z == 0
			α2sπ = v.α2 / √π
			return v.CoulombScaleϵ * NS * (qsz * erfc(qsz2α) - α2sπ * exp(-qsz2α^2) + α2sπ) / π
		else
			eqz = exp(qsz * z)
			αz = v.α * z
			return v.CoulombScaleϵ * NS * (erfc(qsz2α + αz) * eqz - erfc(qsz2α - αz) / eqz + 2 * erf(αz)) / (2π * z)
		end
	else
		error("Wrong method for a ReciprocalGauss3D.")
	end
end
function ReciprocalGauss2D(; ϵ::Real = 1.0, α::Real = 1.0, S::Real = 1.0)
	CoulombScaleϵ = CoulombScale * π / (S * ϵ)
	α2 = 2 * α
	return ReciprocalGauss2D(ϵ, α, S, CoulombScaleϵ, α2)
end

# struct ReciprocalGaussAnisotropy3D <: ReciprocalCoulomb
# 	ϵ::Mat3{Float64}
# 	α::Float64
# 	Ω::Float64
# 	CoulombScale::Float64
# 	α²4::Float64
# end
# function (v::ReciprocalGaussAnisotropy3D)(k::AbstractVector{<:Real})::Float64
# 	k2 = k[1] * k[1] + k[2] * k[2] + k[3] * k[3]
# 	return k2 == 0 ? 0.0 : v.CoulombScale * exp(-k2 / v.α²4) / (transpose(k) * v.ϵ * k)
# end
# function ReciprocalGaussAnisotropy3D(; ϵ::AbstractMatrix{<:Real} = [1 0 0; 0 1 0; 0 0 1], α::Real = 1, Ω::Real = 1)
# 	CoulombScale = CoulombScale * 4π / Ω
# 	α²4 = 4 * α^2
# 	return ReciprocalGaussAnisotropy3D(ϵ, α, Ω, CoulombScale, α²4)
# end
