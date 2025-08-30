export kgridmap
function kgridmap(kgrid::MonkhorstPack, mode = +)

	if !iszero(kgrid.kshift)
		error("Only work for Γ-centered kgrid!")
	end

	redkgrid = RedKgrid(kgrid)
	return _kgrid_map(redkgrid.kdirect, redkgrid.kdirect, mode)
end
function kgridmap(redkgrid::RedKgrid, mode = +)

	if !iszero(irredkgrid.kshift)
		error("Only work for Γ-centered kgrid!")
	end

	return _kgrid_map(redkgrid.kdirect, redkgrid.kdirect, mode)
end

function kgridmap(irredkgrid::IrredKgrid, mode = +)

	if !iszero(irredkgrid.kshift)
		error("Only work for Γ-centered kgrid!")
	end

	return _kgrid_map(irredkgrid.redkdirect, irredkgrid.redkdirect, mode)
end

function kgridmap(aimkdirects, kdirects, mode)

	Nk = length(kdirects)

	mapmatrix = Matrix{Int}(undef, Nk, Nk)
	Threads.@threads for I in CartesianIndices(mapmatrix)
		(i, j) = Tuple(I)
		Tk = mode(kdirects[i], kdirects[j])
		mapmatrix[I] = findfirst(k -> all(isinteger, k - Tk), aimkdirects)
	end

	return mapmatrix
end
