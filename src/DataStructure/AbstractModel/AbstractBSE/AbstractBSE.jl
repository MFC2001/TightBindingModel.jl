
abstract type AbstractBSE <: AbstractModel end

include("./UwithLR.jl")

include("./BSE.jl")
include("./BSEqgrid.jl")

include("./shiftGamma.jl")
include("./elephase.jl")

include("./Kernal.jl")
include("./Hamilton.jl")

function BSEwannier(TB::AbstractTightBindModel, qgrid::MonkhorstPack, paras...; kwards...)

	p = count(TB.period)
	if p == 3
		kgrid = qgrid.kgrid_size .* 2
		kshift = [1 // 2, 1 // 2, 1 // 2]
	elseif p == 2
		kgrid = collect(qgrid.kgrid_size .* 2)
		kshift = [1 // 2, 1 // 2, 1 // 2]
		if !TB.period[1]
			kgrid[1] = 1
			kshift[1] = 0
		elseif !TB.period[2]
			kgrid[2] = 1
			kshift[2] = 0
		elseif !TB.period[3]
			kgrid[3] = 1
			kshift[3] = 0
		end
	elseif p == 1
		if TB.period[1]
			kgrid = [qgrid.kgrid_size[1] * 2, 1, 1]
			kshift = [1 // 2, 0, 0]
		elseif TB.period[2]
			kgrid = [1, qgrid.kgrid_size[1] * 2, 1]
			kshift = [0, 1 // 2, 0]
		elseif TB.period[3]
			kgrid = [1, 1, qgrid.kgrid_size[1] * 2]
			kshift = [0, 0, 1 // 2]
		end
	elseif p == 0
		kgrid = [1, 1, 1]
		kshift = [0, 0, 0]
	end
	#basis is twice as dense
	kgrid = MonkhorstPack(kgrid; kshift)

	return BSEqgrid(TB, kgrid, paras...; kwards...)
end

function (bse::AbstractBSE)(q::ReducedCoordinates)
	if bse.type == :NP
		bandkq, Kᵈ, Kˣ = Kernal!(bse, q)
		return BSE_NP_Hamilton(bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)
	elseif bse.type == :SP
		bandkq, Kᵈ, Kˣ = Kernal!(bse, q)
		return BSE_SP_Hamilton(bse.vckmap, bse.bandk, bandkq, Kᵈ, Kˣ, bse.scissor)
	end
end

# function Base.setproperty!(bse::BSE, prop::Symbol, val)
# 	if prop == :kgrid
# 		if val isa MonkhorstPack
# 			val = RedKgrid(val)
# 		elseif val isa RedKgrid
# 		else
# 			error("Wrong kgrid to BSE.")
# 		end
# 		setfield!(bse, prop, val)
# 		vckmap = vckMap(bse.vckmap.idx2v, bse.vckmap.idx2c, length(bse.kgrid))
# 		setfield!(bse, :vckmap, vckmap)
# 		bandk = BAND(val, bse.TB; vector = true)
# 		setfield!(bse, :bandk, bandk)
# 	elseif prop == :v
# 		vckmap = vckMap(prop, bse.vckmap.idx2c, length(bse.kgrid))
# 		setfield!(bse, :vckmap, vckmap)
# 	elseif prop == :c
# 		vckmap = vckMap(bse.vckmap.idx2v, prop, length(bse.kgrid))
# 		setfield!(bse, :vckmap, vckmap)
# 	elseif prop == :vckmap
# 		if prop.nk ≠ length(bse.kgrid)
# 			error("Please check the vckmap.")
# 		end
# 		setfield!(bse, :vckmap, prop)
# 	elseif prop == :bandk
# 		error("Please don't change BSE.vckmap or BSE.bandk directly.")
# 	else
# 		setfield!(bse, prop, val)
# 	end
# end




