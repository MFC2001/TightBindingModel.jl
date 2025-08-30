
include("./NP.jl")
include("./SP.jl")
include("./wannier_pp.jl")

include("./mmn.jl")
include("./amn.jl")
include("./guess.jl")

export BSE_wannier

function BSE_wannier(qgrid::MonkhorstPack, bse::AbstractBSE; kwards...)
	return BSE_wannier(RedKgrid(qgrid), bse; kwards...)
end
function BSE_wannier(qgrid::RedKgrid, bse::AbstractBSE; kwards...)
	if bse.type == :NP
		return BSE_NP_wannier(qgrid, bse; kwards...)
	elseif bse.type == :SP
		return BSE_SP_wannier(qgrid, bse; kwards...)
	else
		error("Wrong BSE type!")
	end
end
