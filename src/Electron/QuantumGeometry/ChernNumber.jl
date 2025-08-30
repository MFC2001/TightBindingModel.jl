export ChernNumber

function ChernNumber(kgrid, TB::AbstractTightBindModel, n = 1)

	QG = QuantumGeometry(kgrid, TB, n)
	BC = BerryCurve(QG)

	𝐚, 𝐛, 𝐜 = basisvectors(reciprocal(TB.lattice))
	if kgrid.kgrid_size[1] == 1
		S = abs((𝐛×𝐜)[3])
		BC = BC[:, 2, 3, :, :]
	elseif kgrid.kgrid_size[2] == 1
		S = abs((𝐚×𝐜)[3])
		BC = BC[:, 1, 3, :, :]
	elseif kgrid.kgrid_size[3] == 1
		S = abs((𝐚×𝐛)[3])
		BC = BC[:, 1, 2, :, :]
	else
		error("Only calculate Chern Number for 2D.")
	end

	if n isa Integer
		return sum(BC) * S / length(kgrid) / 2π
	elseif n isa AbstractVector
		TBC = similar(BC, size(BC, 1))
		Threads.@threads for k in axes(BC, 1)
			TBC[k] = tr(BC[k, :, :])
		end
		return sum(TBC) * S / length(kgrid) / 2π
	end
end
