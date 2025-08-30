struct BaseReciprocalHoppings{T} <: ReciprocalHoppings{T}
	norb::Int
	Nhop::Matrix{Int}
	hops::Matrix{Vector{Hopping{T}}}
end
function BaseReciprocalHoppings(hr::HR{T}) where {T}
	norb = numorb(hr)
	if sort(hr.orbindex) == 1:norb
		hops = _buildhops_from_hr(norb, hr.path, hr.value)
		Nhop = length.(hops)
	else
		error("Won't construct `ReciprocalHoppings from a nonstandard `HR`.")
	end
	return BaseReciprocalHoppings{T}(norb, Nhop, hops)
end
