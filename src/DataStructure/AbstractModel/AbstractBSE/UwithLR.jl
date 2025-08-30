
#Note ϕK(0) is redefined, is not divergent.
#But we need its divergence, so we can't input k=0.

abstract type GaussLRCorrection end

function (v::GaussLRCorrection)(k::AbstractVector)
	A = zeros(ComplexF64, v.norb, v.norb)
	return v(Val(+), A, k)
end
function (v::GaussLRCorrection)(A, k::AbstractVector)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	A .= 0
	return v(Val(+), A, k)
end
function (v::GaussLRCorrection)(sym::Symbol, nk::Integer)
	if sym ≠ :head
		error("Wrong method for a GaussLRCorrection.")
	end
	A = zeros(ComplexF64, v.norb, v.norb)
	return v(Val(+), A, sym, nk)
end
function (v::GaussLRCorrection)(A, sym::Symbol, nk::Integer)
	if sym ≠ :head
		error("Wrong method for a GaussLRCorrection.")
	end
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	A .= 0
	return v(Val(+), A, sym, nk)
end
function (v::GaussLRCorrection)(k::AbstractVector, nk::Integer)
	A = zeros(ComplexF64, v.norb, v.norb)
	return v(Val(+), A, k, nk)
end
function (v::GaussLRCorrection)(A, k::AbstractVector, nk::Integer)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	A .= 0
	return v(Val(+), A, k, nk)
end

struct GaussLRCorrection3D <: GaussLRCorrection
	φ::ReciprocalGauss3D
	rlattice::ReciprocalLattice
	NG::Int
	G_frac::Vector{Vec3{Int}}
	G_car::Vector{Vec3{Float64}}
	G_orb_phase::Array{ComplexF64}
	norb::Int
	Δorb::Matrix{Vec3{Float64}}
end
function (v::GaussLRCorrection3D)(::Val{+}, A, k::AbstractVector, nk::Integer)
	if iszero(k)
		return v(Val(+), A, :head, nk)
	else
		return v(Val(+), A, k)
	end
end
function (v::GaussLRCorrection3D)(::Val{+}, A, k::AbstractVector)

	k_car = v.rlattice * k
	φkG = map(G -> v.φ(k_car + G), v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = cis(-2π * (k ⋅ v.Δorb[i, j])) * sum(Gi -> φkG[Gi] * v.G_orb_phase[Gi, i, j], Base.OneTo(v.NG))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	Vᵢᵢ = sum(φkG)
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function (v::GaussLRCorrection3D)(::Val{+}, A, sym::Symbol, nk)

	G0 = findfirst(iszero, v.G_frac)
	I = setdiff(1:v.NG, G0)

	φG = map(G -> v.φ(G), v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = sum(Gi -> φG[Gi] * v.G_orb_phase[Gi, i, j], I)
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	Vᵢᵢ = sum(φG)
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	φhead = v.φ(sym, nk)
	A .+= φhead

	return A
end
function GaussLRCorrection3D(lattice::Lattice, orblocat_frac::AbstractVector, α::Real; δ = 1e-6, ϵ = 1)

	G2_max = -4 * α^2 * log(δ) * 1.2

	rlattice = reciprocal(lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = ((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)

	Ggrid = Int.(cld.(√G2_max * 1.1, [h₁, h₂, h₃])) * 2 .+ 1
	Ggrid = gridindex(Ggrid)
	G_frac = filter(G -> sum(abs2, rlattice * G) < G2_max, Ggrid)
	G_car = map(G -> rlattice * G, G_frac)

	Ω = abs((lattice[:, 1] × lattice[:, 2]) ⋅ lattice[:, 3])
	φ = ReciprocalGauss3D(; ϵ, α, Ω)

	norb = length(orblocat_frac)
	Δorb = Matrix{Vec3{Float64}}(undef, norb, norb)
	NG = length(G_frac)
	G_orb_phase = Array{ComplexF64}(undef, NG, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocat_frac[j] - orblocat_frac[i]
		Δorb[i, j] = dorb
		Δorb[j, i] = -dorb
		for (Gi, G) in enumerate(G_frac)
			G_orb_phase[Gi, i, j] = cis(-2π * (G ⋅ Δorb[i, j]))
		end
		G_orb_phase[:, j, i] = conj.(G_orb_phase[:, i, j])
	end

	return GaussLRCorrection3D(φ, rlattice, NG, G_frac, G_car, G_orb_phase, norb, Δorb)
end

"""
	输入的k可能会在非周期方向非0，例如计算激子wannier时，但按照Poisson求和得到的公式中并不包含非周期方向的值
"""
struct GaussLRCorrection2D <: GaussLRCorrection
	φ::ReciprocalGauss2D
	rlattice::ReciprocalLattice
	NG::Int
	G_frac::Vector{Vec3{Int}}
	G_car::Vector{Vec3{Float64}}
	G_orb_phase::Array{ComplexF64}
	norb::Int
	Δorb::Matrix{Vec3{Float64}}
	Δorbz::Matrix{Float64}
	xyindex::SVector{2, Int}
	zindex::Int
end
function (v::GaussLRCorrection2D)(::Val{+}, A, k::AbstractVector, nk::Integer)
	#Only use kx,ky.
	k2D = k[v.xyindex]

	if iszero(k2D)
		return v(Val(+), A, :head, nk)
	else
		return v(Val(+), A, k)
	end
end
function (v::GaussLRCorrection2D)(::Val{+}, A, k::AbstractVector)

	#Only use kx,ky.
	k2D = [0.0, 0.0, 0.0]
	k2D[v.xyindex] = k[v.xyindex]

	k_car = v.rlattice * k2D
	kG_car = map(G -> k_car + G, v.G_car)

	for j in 2:v.norb, i in 1:j-1
		vlr = cis(-2π * (k2D ⋅ v.Δorb[i, j])) * sum(Gi -> v.φ(kG_car[Gi], v.Δorbz[i, j]) * v.G_orb_phase[Gi, i, j], Base.OneTo(v.NG))
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	Vᵢᵢ = sum(Gi -> v.φ(kG_car[Gi], 0), Base.OneTo(v.NG))
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function (v::GaussLRCorrection2D)(::Val{+}, A, sym::Symbol, nk::Integer)

	G0 = findfirst(iszero, v.G_frac)
	I = setdiff(1:v.NG, G0)

	for j in 2:v.norb, i in 1:j-1
		φhead = v.φ(sym, nk, v.Δorbz[i, j])
		vlr = sum(Gi -> v.φ(v.G_car[Gi], v.Δorbz[i, j]) * v.G_orb_phase[Gi, i, j], I) + φhead
		A[i, j] += vlr
		A[j, i] += conj(vlr)
	end
	φhead = v.φ(sym, nk, 0)
	Vᵢᵢ = sum(Gi -> v.φ(v.G_car[Gi], 0), I) + φhead
	for i in 1:v.norb
		A[i, i] += Vᵢᵢ
	end

	return A
end
function GaussLRCorrection2D(lattice::Lattice, orblocat_frac::AbstractVector, period::AbstractVector, α::Real; δ = 1e-6, ϵ = 1)

	G2_max = -4 * α^2 * log(δ) * 1.5

	rlattice = reciprocal(lattice)
	b₁ = rlattice[:, 1]
	b₂ = rlattice[:, 2]
	b₃ = rlattice[:, 3]

	V_BZ = ((b₁ × b₂) ⋅ b₃)

	h₁ = V_BZ / norm(b₂ × b₃)
	h₂ = V_BZ / norm(b₃ × b₁)
	h₃ = V_BZ / norm(b₁ × b₂)

	Ggrid = Int.(cld.(√G2_max * 1.2, [h₁, h₂, h₃])) * 2 .+ 1

	a₁ = lattice[:, 1]
	a₂ = lattice[:, 2]
	a₃ = lattice[:, 3]

	if !period[1]
		Ggrid[1] = 1
		S = norm(a₂ × a₃)
		xyindex = SVector{2, Int}(2, 3)
		zindex = 1
	elseif !period[2]
		Ggrid[2] = 1
		S = norm(a₃ × a₁)
		xyindex = SVector{2, Int}(1, 3)
		zindex = 2
	elseif !period[3]
		Ggrid[3] = 1
		S = norm(a₁ × a₂)
		xyindex = SVector{2, Int}(1, 2)
		zindex = 3
	end

	Ggrid = gridindex(Ggrid)
	G_frac = filter(G -> sum(abs2, rlattice * G) < G2_max, Ggrid)
	G_car = map(G -> rlattice * G, G_frac)

	φ = ReciprocalGauss2D(; ϵ, α, S)

	norb = length(orblocat_frac)
	Δorb = Matrix{Vec3{Float64}}(undef, norb, norb)
	NG = length(G_frac)
	G_orb_phase = Array{ComplexF64}(undef, NG, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocat_frac[j] - orblocat_frac[i]
		Δorb[i, j] = dorb
		Δorb[j, i] = -dorb
		for (Gi, G) in enumerate(G_frac)
			G_orb_phase[Gi, i, j] = cis(-2π * (G ⋅ Δorb[i, j]))
		end
		G_orb_phase[:, j, i] = conj.(G_orb_phase[:, i, j])
	end

	Δorbz = map(dorb -> (lattice*dorb)[zindex], Δorb)

	return GaussLRCorrection2D(φ, rlattice, NG, G_frac, G_car, G_orb_phase, norb, Δorb, Δorbz, xyindex, zindex)
end

#TODO
struct GaussLRCorrection1D <: GaussLRCorrection
	norb::Int
end
function (v::GaussLRCorrection1D)(::Val{+}, A, para...)
	error("TODO")
	return A
end
struct GaussLRCorrection0D <: GaussLRCorrection
	norb::Int
end
function (v::GaussLRCorrection0D)(::Val{+}, A, para...)
	return A
end


"""
	Direct coulomb term with long-range correction from Gauss potential in reciprocal space.
"""
# We calculate head term instead of k=0.
# 0D情况，倘若U是实的，那么计算结果也是实的，其余情况均为复的，这种情况特殊处理，均归于BSEcluster
# 这里的0D长程修正仅用于不考虑长程修正的情况。
struct UwithLR{T <: Number, V <: GaussLRCorrection}
	norb::Int
	SR::BaseReciprocalHoppings{T}
	LR::V
end
function (v::UwithLR)(parameters...)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, parameters...)
end
function (v::UwithLR)(A::AbstractMatrix, k::AbstractVector, nk::Integer)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.SR(A, k)
	v.LR(Val(+), A, k, nk)
	return A
end
function (v::UwithLR)(A::AbstractMatrix, k::AbstractVector)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.SR(A, k)
	v.LR(Val(+), A, k)
	return A
end
function (v::UwithLR)(A::AbstractMatrix, sym::Symbol, nk::Integer)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.SR(A, [0, 0, 0])
	v.LR(Val(+), A, sym, nk)
	return A
end

function UwithLR(TB, U, α, δ, ϵ)

	p = count(TB.period)

	if p == 0
		error("Please use BSEcluster")
	end

	V_SR = _UwithLR_SR(TB, U, α, ϵ)
	if p == 3
		V_LR = GaussLRCorrection3D(TB.lattice, TB.orb_location, α; δ, ϵ)
	elseif p == 2
		V_LR = GaussLRCorrection2D(TB.lattice, TB.orb_location, TB.period, α; δ, ϵ)
	elseif p == 1
		#TODO 
		error("TODO")
		# V_LR = GaussLRCorrection1D(TB.lattice, TB.orb_location, TB.period, α; δ, ϵ)
	end

	return UwithLR(numorb(TB), V_SR, V_LR)
end

function _UwithLR_SR(TB, U, α, ϵ)

	φR = RealGauss(; ϵ, α)

	hops_SR = deepcopy(U.hops)

	for I in CartesianIndices(U.hops)
		(i, j) = Tuple(I)
		for (ii, hop) in enumerate(U.hops[I])
			value_SR = hop.t - φR(TB.lattice * (hop.R + TB.orb_location[j] - TB.orb_location[i]))
			hops_SR[I][ii] = similar(hop, value_SR)
		end
	end

	return BaseReciprocalHoppings(U.norb, U.Nhop, hops_SR)
end


function _UwithLR_0D(TB, U, ϵ)

	norb = numorb(U)
	hops = _buildhops_RH(norb, U.path, U.value)
	Nhop = length.(hops)

	I = findfirst(n -> n > 1, Nhop)
	if !isnothing(I)
		error("Wrong 0D U_hr.")
	end

	VR = V3DR(; ϵ)

	I = findall(iszero, Nhop)
	for II in I
		(i, j) = Tuple(II)
		push!(hops[II], Hopping{Float64}(i, j, [0, 0, 0], VR(TB.lattice * (TB.orb_location[j] - TB.orb_location[i]))))
		Nhop[II] = 1
	end

	SR = BaseReciprocalHoppings(norb, hops, Nhop)
	LR = GaussLRCorrection0D(norb)

	return UwithLR(norb, SR, LR)
end

function UwithLR(SR::ReciprocalHoppings{T}) where {T}
	#This means you don't want include long-range correction.
	SR = BaseReciprocalHoppings(SR.norb, SR.hops, SR.Nhop)
	LR = GaussLRCorrection0D(SR.norb)
	return UwithLR{T}(SR.norb, SR, LR)
end
