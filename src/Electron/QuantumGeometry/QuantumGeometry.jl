
include("./abelian.jl")
include("./nonabelian.jl")
include("./partialHamilton.jl")
include("./ChernNumber.jl")

export QuantumGeometry, BerryCurve, QuantumMetric

function QuantumGeometry(kgrid, TB::AbstractTightBindModel, n = 1)
	# n is bandindex
	hr = HR(eachindex(TB.orb_location), TB.hop, TB.Nhop, Matrix{T}(undef, 0, 0), Vector{T}(undef, 0))

	if kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid).kdirect
	elseif kgrid isa RedKgrid || kgrid isa IrredKgrid
		kgrid = kgrid.kdirect
	end

	return QuantumGeometry(kgrid, hr, TB.lattice, TB.orb_location, n)
end

function BerryCurve(QG::AbstractArray)
	if ndims(QG) == 3
		return imag.(conj.(QG) .- QG)
	elseif ndims(QG) == 5
		#F is the antisymmetric part of QG.
		F = similar(QG, ComplexF64)
		for k in axes(QG, 1), i in 1:3, j in 1:3
			T = QG[k, i, j, :, :]
			F[k, i, j, :, :] = 1im * (T .- T')
		end
		return F
	end
end

function QuantumMetric(QG::AbstractArray)
	if ndims(QG) == 3
		return real.(conj.(QG) .+ QG) ./ 2
	elseif ndims(QG) == 5
		#g is the symmetric part of QG
		g = similar(QG, ComplexF64)
		for k in axes(QG, 1), i in 1:3, j in 1:3
			T = QG[k, i, j, :, :]
			g[k, i, j, :, :] = (T .+ T') ./ 2
		end
		return g
	end
end
