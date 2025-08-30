
export AbstractModel, TightBindModel, SymTightBindModel, BSE

abstract type AbstractModel end


include("./AbstractTightBindModel.jl")
include("./AbstractBSE/AbstractBSE.jl")
# include("./BSE_supercell.jl")



# function Base.filter(f, symhr::SymHR{T}) where {T}

# 	symhr = deepcopy(symhr)
# 	for I in CartesianIndices(symhr.hop)
# 		symhr.hop[I] = filter(hop -> f([hop.R..., hop.i, hop.j], symhr.value[hop.t]), symhr.hop[I])
# 	end

# 	return symhr
# end

# struct MeanField{UT <: Number, HT <: Number} <: AbstractModel
# 	irredkgrid::IrredKgrid{UT}
# 	symTB::SymTightBindModel{HT}
# 	DensityMatrix::Function
# 	Vik::Function
# end

