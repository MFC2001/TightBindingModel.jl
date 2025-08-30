
abstract type AbstractReciprocalHoppings end

include("./ReciprocalHoppings.jl")
include("./BaseReciprocalHoppings.jl")
include("./HermitianReciprocalHoppings.jl")

# include("./GammaReciprocalHoppings.jl")




# struct ReciprocalHoppingsAtomic{T <: Number} <: AbstractReciprocalHoppings{T}
# 	hops::Matrix{Vector{Hopping{T}}}
# 	Nhop::Matrix{Int}
# 	orblocat::Vector{Vec3{Float64}}
# end
# function ReciprocalHoppingsAtomic(hr::HR{T}, orblocat::AbstractVector) where {T}
# 	norb = numorb(hr)
# 	if sort(hr.orbindex) == 1:norb
# 		hops, Nhop = _build_RH(norb, path, value)
# 	else
# 		error("Won't construct `ReciprocalHoppings from a nonstandard `HR`.")
# 	end
# 	return ReciprocalHoppingsAtomic{T}(hops, Nhop, orblocat)
# end

# @inline numorb(rh::ReciprocalHoppingsAtomic{T}) = size(rh.hops, 1) where {T}
# @inline function (rh::ReciprocalHoppingsAtomic{T})(k) where {T}
# 	norb = numorb(rh)
# 	A = Matrix{ComplexF64}(undef, norb, norb)
# 	return rh(A, k)
# end
# @inline function (rh::ReciprocalHoppingsAtomic{T})(A, k) where {T}
# 	size(A) == size(rh.hops) || error("Buffer size mismatch.")
# 	map(CartesianIndices(rh.hops)) do I
# 		A[I] = iszero(rh.Nhop[I]) ? 0 : cis(2π * (k ⋅ (rh.orblocat[j] - rh.orblocat[i]))) * sum(hop -> hop.t * cis(2π * (k ⋅ hop.R)), rh.hops[I])
# 	end
# 	return A
# end


