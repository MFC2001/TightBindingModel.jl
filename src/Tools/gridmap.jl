function gridmap(grid::AbstractVector, mode = +)

	N = length(grid)

	N1_max = maximum(x -> x[1], grid)
	N1_min = minimum(x -> x[1], grid)
	N2_max = maximum(x -> x[2], grid)
	N2_min = minimum(x -> x[2], grid)
	N3_max = maximum(x -> x[3], grid)
	N3_min = minimum(x -> x[3], grid)

	N1 = N1_max - N1_min + 1
	N2 = N2_max - N2_min + 1
	N3 = N3_max - N3_min + 1

	Nv = Vec3(N1, N2, N3)

	frac_grid = map(x -> Vec3(x .// Nv), grid)

	mapmatrix = Matrix{Int}(undef, N, N)
	Threads.@threads for I in CartesianIndices(mapmatrix)
		(i, j) = Tuple(I)
		T = mode(frac_grid[i], frac_grid[j])
		mapmatrix[I] = findfirst(x -> all(isinteger, x - T), frac_grid)
	end

	return mapmatrix
end
