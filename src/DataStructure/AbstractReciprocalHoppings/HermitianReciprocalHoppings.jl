
struct HermitianReciprocalHoppings{T <: Number, U <: ReciprocalHoppings{T}} <: AbstractReciprocalHoppings
	rh::U
	uplo::Symbol
end
# const SingleParticalHamiltonian = HermitianReciprocalHoppings
@inline function (hrh::HermitianReciprocalHoppings{T, U})(parameters...) where {T <: Number, U <: ReciprocalHoppings{T}}
	A = Matrix{ComplexF64}(undef, hrh.rh.norb, hrh.rh.norb)
	return hrh(A, parameters...)
end
@inline function (hrh::HermitianReciprocalHoppings)(A::AbstractMatrix, parameters...)
	for i in Base.OneTo(hrh.rh.norb)
		A[i, i] = hrh.rh(i, i, parameters...)
	end
	_hrh_offdiag_term!(Val(hrh.uplo), A, hrh.rh, parameters...)
	return Hermitian(A, hrh.uplo)
end
@inline function _hrh_offdiag_term!(::Val{:U}, A, rh, parameters...)
	for j in 2:rh.norb, i in 1:j-1
		A[i, j] = rh(i, j, parameters...)
	end
end
@inline function _hrh_offdiag_term!(::Val{:L}, A, rh, parameters...)
	for j in 1:rh.norb, i in j+1:rh.norb
		A[i, j] = rh(i, j, parameters...)
	end
end

HR(hrh::HermitianReciprocalHoppings) = HR(hrh.rh)
numorb(hrh::HermitianReciprocalHoppings) = hrh.rh.norb

function Base.show(io::IO, hrh::HermitianReciprocalHoppings)
	print(io, "HermitianReciprocalHoppings with $(sum(hrh.rh.Nhop)) hoppings and $(hrh.rh.norb) orbitals.")
end

function HermitianReciprocalHoppings(hr::HR{T}, uplo::Symbol = :U) where {T}
	if uplo ∉ (:U, :L)
		throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
	end
	rh = BaseReciprocalHoppings(hr)
	return HermitianReciprocalHoppings(rh, uplo)
end
function LinearAlgebra.Hermitian(rh::U, uplo::Symbol = :U) where {T <: Number, U <: ReciprocalHoppings{T}}
	if uplo ∉ (:U, :L)
		throw(ArgumentError("uplo must be either :U or :L, got :$uplo"))
	end
	return HermitianReciprocalHoppings{T, U}(rh, uplo)
end

Base.zero(::HermitianReciprocalHoppings{T, U}, norb::Integer) where {T <: Number, U <: ReciprocalHoppings{T}} =
	Hermitian(zero(U, norb), :U)
Base.zero(::Type{HermitianReciprocalHoppings{T, U}}, norb::Integer) where {T <: Number, U <: ReciprocalHoppings{T}} =
	Hermitian(zero(U, norb), :U)
