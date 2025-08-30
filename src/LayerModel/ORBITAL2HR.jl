export ORBITAL2HR
function ORBITAL2HR(poscar::POSCAR, orbital₁::ORBITAL, orbindex₁::AbstractVector{<:Integer}, orbital₂::ORBITAL, orbindex₂::AbstractVector{<:Integer}, V::Function; heps::Real = 1e-6, maxdist = findmaxdistV(V, heps))::HR

	#Expand supercell several times in different cases.
	times = POSCAR2HRs.expandtimes(poscar.lattvec, maxdist, poscar.periodicity)
	supercell₂ = POSCAR2HRs.SuperCell(orbital₂.location, poscar.lattvec, poscar.periodicity; atom_index = orbindex₂)
	POSCAR2HRs.ExpandSuperCell!(supercell₂, times)

	heps2 = heps^2
	maxdist2 = maxdist^2

	hr_path = Vector{SVector{3, Int}}(undef, 0)
	hr_value = Vector{Float64}(undef, 0)
	hr_orbpath = Matrix{Int}(undef, 2, 0)

	lk = ReentrantLock()
	Threads.@threads for i in eachindex(orbital₁.location)
		startindex = orbindex₁[i]
		D = supercell₂.location .- [orbital₁.location[i]]
		d = map(norm2, D)
		# I = findall(d .< maxdist2)

		I = d .< maxdist2
		value = V.(D[I])
		path = supercell₂.path[I]
		orbpath = transpose([fill(startindex, length(value)) supercell₂.index[I]])

		II = abs2.(value) .> heps2
		value = value[II]
		path = path[II]
		orbpath = orbpath[:, II]

		value = [value; conj.(value)]
		path = [path; -path]
		orbpath = [orbpath orbpath[[2, 1], :]]

		# @show I
		lock(lk) do
			append!(hr_path, path)
			append!(hr_value, value)
			hr_orbpath = [hr_orbpath orbpath]
		end
	end

	return HR(hr_path, hr_orbpath, hr_value; buildindex = "N")
end
