function atompath_sum(location, supercell, atompath::AbstractVector{AtomPath}, unitlattice, maxlattdist::Real, deps::Real)
	pathvec = Vector{SVector{3, Float64}}(undef, length(atompath))
	for i in eachindex(atompath)
		pathvec[i] = unitlattice * atompath[i].unitpath
	end

	maxlattdist2 = maxlattdist^2
	deps2 = deps^2

	atompaths = Vector{AtomPath}(undef, 0)
	atompathindex = Vector{SVector{2, Int}}(undef, 0)


	for i in eachindex(location)

		D = map(x -> x - location[i], supercell.location)
		d = map(norm2, D)
		#P.S. location has contained the offset between atoms.
		I = findall(d .< maxlattdist2)

		T = map(x -> findmin(t -> norm2(t - x), pathvec), D[I])
		if any(t -> t[1] > deps2, T)
			# maximum([t[1] for t in T]) > deps2
			error("Wrong poscar's path from atompath_sum in poscar2hr.")
		end

		T = [(I[i], T[i][2]) for i in eachindex(I)]
		atompaths_add = [AtomPath(supercell.path[II], atompath[aim].hrindex) for (II, aim) in T]
		atompathindex_add = [SVector{2}(i, supercell.index[II]) for (II, ~) in T]

		append!(atompaths, atompaths_add)
		append!(atompathindex, atompathindex_add)

		# for II in I
		# 	(min_d, aim) = findmin(norm2, pathvec .- [D[II]])
		# 	if min_d > deps2
		# 		error("Wrong poscar's path from atompath_sum in poscar2hr.")
		# 	end
		# 	push!(atompaths, AtomPath(supercell.path[II], atompath[aim].hrindex))
		# 	push!(atompathindex, SVector{2}(i, supercell.index[II]))
		# end
	end

	N_atompaths = map(ap -> length(ap.hrindex), atompaths)
	hr_path_unithrindex = Matrix{Int}(undef, 6, sum(N_atompaths))
	n = 0
	for (i, ap) in enumerate(atompaths)
		if N_atompaths[i] > 0
			hr_path_unithrindex[:, n+1:n+N_atompaths[i]] = [repeat([ap.unitpath; atompathindex[i]], inner = (1, N_atompaths[i])); transpose(ap.hrindex)]
			n += N_atompaths[i]
		end
	end




	# hr_path_sum = Vector{SVector{3, Int}}(undef, 0)
	# hr_orbpath_sum = Matrix{Int}(undef, 3, 0)
	# for i in eachindex(atompaths)
	# 	ap = atompaths[i]
	# 	n = length(ap.hrindex)
	# 	if n > 0
	# 		append!(hr_path_sum, fill(ap.unitpath, n))
	# 		hr_orbpath_sum = [hr_orbpath_sum [repeat(atompathindex[i], inner = (1, n)); transpose(ap.hrindex)]]
	# 	end
	# end

	return hr_path_unithrindex
end

function atompath_sum_threads(location, supercell, atompath::AbstractVector{AtomPath}, unitlattice, maxlattdist::Real, deps::Real)
	pathvec = Vector{SVector{3, Float64}}(undef, length(atompath))
	for i in eachindex(atompath)
		pathvec[i] = unitlattice * atompath[i].unitpath
	end

	maxlattdist = maxlattdist^2
	deps2 = deps^2

	atompaths = Vector{AtomPath}(undef, 0)
	atompathindex = Vector{SVector{2, Int}}(undef, 0)

	lk = ReentrantLock()
	Threads.@threads for i in eachindex(location)

		D = map(x -> x - location[i], supercell.location)
		d = map(norm2, D)
		#P.S. location has contained the offset between atoms.
		I = findall(d .< maxlattdist)

		T = map(x -> findmin(norm2, pathvec .- [x]), D[I])
		if maximum([t[1] for t in T]) > deps2
			error("Wrong poscar's path from atompath_sum in poscar2hr.")
		end

		T = [(I[i], T[i][2]) for i in eachindex(I)]
		atompaths_add = [AtomPath(supercell.path[II], atompath[aim].hrindex) for (II, aim) in T]
		atompathindex_add = [SVector{2}(i, supercell.index[II]) for (II, ~) in T]

		lock(lk) do
			append!(atompaths, atompaths_add)
			append!(atompathindex, atompathindex_add)
		end
	end

	hr_path_sum = Vector{SVector{3, Int}}(undef, 0)
	hr_orbpath_sum = Matrix{Int}(undef, 3, 0)
	for i in eachindex(atompaths)
		ap = atompaths[i]
		n = length(ap.hrindex)
		if n > 0
			append!(hr_path_sum, fill(ap.unitpath, n))
			hr_orbpath_sum = [hr_orbpath_sum [repeat(atompathindex[i], inner = (1, n)); transpose(ap.hrindex)]]
		end
	end

	return hr_path_sum, hr_orbpath_sum
end
