export BAND

include("./TightBindingModel/TightBindingModel.jl")
include("./BSE/BSE.jl")


function BAND(kgrid::MonkhorstPack, args...; vector::Bool = false)
	return BAND(RedKgrid(kgrid), args...; vector)
end
function BAND(kpoints::AbstractBrillouinZone, args...; vector::Bool = false)
	return BAND(kpoints.kdirect, args...; vector)
end
