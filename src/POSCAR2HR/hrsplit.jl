function hrsplit(hr::HR, cell::Cell, orbital_belonging::AbstractVector{<:Integer}, deps::Real)::AtomHR
	#Split unit hr to each atom pair's hrs.
	norb = length(orbital_belonging)
	atom_name = cell.name
	atom_num = length(cell.location)

	#Calculate the max distance, and get the logical transitions with distances shorter than maxlattdist.
	(maxlattdist, I) = maxhrpath(hr.path, cell, deps)

	hr_index = [Vector{typeof(length(hr.value))}(undef, 0) for _ in CartesianIndices((norb, norb))]
	for i in axes(hr.path, 1)
		push!(hr_index[hr.path[i, 4], hr.path[i, 5]], i)
	end

	#orbital_index = [[1,2,3],[4,5,6]]...
	orbital_index = [Vector{typeof(length(orbital_belonging))}(undef, 0) for _ in 1:atom_num]
	for i in eachindex(orbital_belonging)
		push!(orbital_index[orbital_belonging[i]], i)
	end


	parindex = Matrix{Tuple{Int, Int}}(undef, atom_num, atom_num)
	parname = Matrix{String}(undef, atom_num, atom_num)
	atompath = Matrix{Vector{AtomPath}}(undef, atom_num, atom_num)
	Natompath = Matrix{typeof(length(hr.value))}(undef, atom_num, atom_num)


	Threads.@threads for (i, j) in [(i, j) for i in 1:atom_num, j in 1:atom_num]
		parindex[i, j] = (i, j)
		parname[i, j] = atom_name[i] * "->" * atom_name[j]

		Tatompath = index2atom(orbital_index[i], orbital_index[j], hr.path, hr_index)
		Natompath[i, j] = length(Tatompath)
		if Natompath[i, j] == 0
			continue
		end

		allpath = [t.unitpath for t in Tatompath]
		II = setdiff(I, allpath)
		for n in eachindex(II)
			push!(Tatompath, AtomPath(II[n], Int[]))
		end
		atompath[i, j] = Tatompath
	end


	hrorbindex = transpose(hr_orbital_reindex(hr.path, orbital_index))

	return AtomHR(parindex, parname, atompath, Natompath, maxlattdist, hrorbindex, orbital_index, length.(orbital_index))
end

function maxhrpath(hr_path, cell::Cell, deps::Real)
	#Calculate the max distance, and get the logical positive lattice vectors with distances shorter than maxlattdist.
	lattice = cell.lattice
	period = cell.period

	allpath = unique(map(x -> SVector{3}(x[1:3]), eachrow(hr_path)))

	p = count(period)
	if p == 3

		maxlattdist = sqrt(maximum(map(x -> norm2(lattice * x), allpath))) + deps

		a₁ = lattice[:, 1]
		a₂ = lattice[:, 2]
		a₃ = lattice[:, 3]

		θ₁ = abs(calc_vecangle(a₁, a₂))
		θ₂ = abs(calc_vecangle(a₂, a₃))
		θ₃ = abs(calc_vecangle(a₃, a₁))

		dV = abs((a₁ × a₂) ⋅ a₃) / (norm(a₁) * norm(a₂) * norm(a₃))
		#It's related to the radius of biggest sphere inside parallelepiped.
		d = (2 * maxlattdist + 0.1) * maximum(sin.([θ₁, θ₂, θ₃])) / dV

		#The lattice points in a parallelepiped.
		n₁ = Int(cld(d / 2, sqrt(sum(abs2, a₁))))
		n₂ = Int(cld(d / 2, sqrt(sum(abs2, a₂))))
		n₃ = Int(cld(d / 2, sqrt(sum(abs2, a₃))))
		I = [repeat(-n₁:n₁, inner = 2 * n₂ + 1) repeat(-n₂:n₂, outer = 2 * n₁ + 1)]
		I = [repeat(I, outer = (2 * n₃ + 1, 1)) repeat(-n₃:n₃, inner = (2 * n₁ + 1) * (2 * n₂ + 1))]
		I = [SVector{3}(I[i, :]) for i in axes(I, 1)]

		D = map(x -> norm2(lattice * x), I)
		D = D .< maxlattdist^2
		I = I[D]

	elseif p == 2
		if !period[1]
			T = 2:3
		elseif !period[2]
			T = [3, 1]
		elseif !period[3]
			T = 1:2
		end

		lattvec1 = MMatrix{3, 3}(zeros(3, 3))
		lattvec1[:, T] .= lattice[:, T]

		maxlattdist = sqrt(maximum(map(x -> norm2(lattvec1 * x), allpath))) + deps

		a₁ = lattvec1[:, T[1]]
		a₂ = lattvec1[:, T[2]]

		θ = abs(calc_vecangle(a₁, a₂))
		#It's related to the radius of tangent circle inside parallelogram.
		d = (2 * maxlattdist + 0.1) / sin(θ)

		#The lattice points in a parallelogram.
		n₁ = Int(cld(d / 2, norm(a₁)))
		n₂ = Int(cld(d / 2, norm(a₂)))
		II = [repeat(-n₁:n₁, inner = 2 * n₂ + 1) repeat(-n₂:n₂, outer = 2 * n₁ + 1)]

		I = Vector{SVector{3, eltype(II)}}(undef, size(II, 1))
		TI = zeros(eltype(II), 3)
		for i in eachindex(I)
			TI[T] .= II[i, :]
			I[i] = SVector{3}(TI)
		end

		D = map(x -> norm2(lattvec1 * x), I)
		D = D .< maxlattdist^2
		I = I[D]

	elseif p == 1
		if period[1]
			T = 1
		elseif period[2]
			T = 2
		elseif period[3]
			T = 3
		end

		lattvec1 = MMatrix{3, 3}(zeros(3, 3))
		lattvec1[:, T] .= lattice[:, T]

		maxlattdist = sqrt(maximum(map(x -> norm2(lattvec1 * x), allpath))) + deps

		a₁ = lattvec1[:, T]


		#It's related to the radius of tangent circle inside parallelogram.
		d = 2 * maxlattdist + 0.1

		#The lattice points in a parallelogram.
		n₁ = Int(cld(d / 2, norm(a₁)))
		II = collect(-n₁:n₁)

		I = Vector{SVector{3, eltype(II)}}(undef, length(II))
		TI = zeros(eltype(II), 3)
		for i in eachindex(I)
			TI[T] = II[i]
			I[i] = SVector{3}(TI)
		end

		D = map(x -> norm2(lattvec1 * x), I)
		D = D .< maxlattdist^2
		I = I[D]

	elseif p == 0
		error("Wrong periodicity of unitcell, it should be periodic in at least one dimension.")
	end

	return maxlattdist - deps / 2, I
end

function index2atom(index4, index5, hr_path, hr_index)::Vector{AtomPath}
	if isempty(index4) || isempty(index5)
		return Vector{AtomPath}(undef, 0)
	end

	aim_hr_index = Vector{typeof(size(hr_path, 1))}(undef, 0)
	for i in index4, j in index5
		append!(aim_hr_index, hr_index[i, j])
	end

	atompath = Vector{AtomPath}(undef, 0)

	aim_hr_path = map(i -> SVector{3}(hr_path[i, 1:3]), aim_hr_index)
	while length(aim_hr_index) > 0
		path = aim_hr_path[1]
		I = findall(x -> x == path, aim_hr_path)
		push!(atompath, AtomPath(path, aim_hr_index[I]))

		resindex = setdiff(eachindex(aim_hr_index), I)
		aim_hr_path = aim_hr_path[resindex]
		aim_hr_index = aim_hr_index[resindex]
	end

	return atompath
end

function hr_orbital_reindex(hr_path, orbital_index)

	hr_orbpath = hr_path[:, 4:5]
	hrorbindex = similar(hr_orbpath)

	for index in orbital_index
		for i in eachindex(index)
			T = hr_orbpath .== index[i]
			hrorbindex[T] .= i
		end
	end
	return hrorbindex
end
