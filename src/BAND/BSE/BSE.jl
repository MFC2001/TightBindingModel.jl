
function _BSE_eigen_init()

	# 检查是否使用KrylovKit
	if isdefined(Main, :KrylovKit)
		include("./core_KrylovKit.jl")
	else
		include("./core_LinearAlgebra.jl")
	end

end
_BSE_eigen_init()

function BAND(qpoints::AbstractVector{<:ReducedCoordinates}, bse::AbstractBSE; vector::Bool = false)
	if bse.type == :NP
		return BSE_NP(qpoints, bse, Val(vector))
	elseif bse.type == :SP
		return BSE_SP(qpoints, bse, Val(vector))
	else
		error("Wrong BSE type!")
	end
end
