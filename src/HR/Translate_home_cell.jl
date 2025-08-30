export Translate_home_cell
"""
	Translate_home_cell(hr::HR, translatepair::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

	translatepair is index => R, the orbital index in zero cell is acctually in R cell. 
"""
function Translate_home_cell(hr::HR, translatepair::Pair{<:Integer, <:AbstractVector{<:Integer}}...)

	index = [Vector{Int}(undef, 0) for _ in CartesianIndices(hr.hop)]
	for i in axes(hr.path, 1)
		push!(index[hr.path[i, 4], hr.path[i, 5]], i)
	end
	Nindex = length.(index)


	hr_path = deepcopy(hr.path)
	hr_value = deepcopy(hr.value)
	norb = numorb(hr)

	for p in translatepair

		for i in setdiff(1:norb, p.first)
			if Nindex[i, p.first] > 0
				for ii in index[i, p.first]
					hr_path[ii, 1:3] = hr_path[ii, 1:3] + p.second
				end
			end
		end
		for i in setdiff(1:norb, p.first)
			if Nindex[p.first, i] > 0
				for ii in index[p.first, i]
					hr_path[ii, 1:3] = hr_path[ii, 1:3] - p.second
				end
			end
		end

	end

	return HR(hr_path, hr_value; hrsort = 'Y')
end
function Translate_home_cell(orbital::ORBITAL, lattice::Lattice)

	orbital = deepcopy(orbital)

	orblocation_frac = map(x -> lattice \ x, orbital.location)
	frac_incell = map(x -> mod.(x, 1), orblocation_frac)
	R = map(i -> round.(Int, orblocation_frac[i] - frac_incell[i]), eachindex(frac_incell))

	translatepair = Vector{Pair{Int, Vec3{Int}}}(undef, 0)
	for i in eachindex(R)
		if !iszero(R[i])
			push!(translatepair, i => Vec3{Int}(R[i]))
			orbital.location[i] = lattice * (orblocation_frac[i] - R[i])
		end
	end

	return orbital, translatepair
end
