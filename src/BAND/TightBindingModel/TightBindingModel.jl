
include("./core.jl")

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
