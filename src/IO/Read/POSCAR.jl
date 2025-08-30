export ReadPOSCAR

"""
	ReadPOSCAR(file::AbstractString; period = Bool[1, 1, 1])::Cell

Read POSCAR, return a Cell.
"""
function ReadPOSCAR(file::AbstractString; period = [1, 1, 1])::Cell

	file = open(file, "r")
	readline(file)
	scale = parse(Float64, readline(file))
	lattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		lattice[i, :] = [parse(Float64, ss) for ss in split(readline(file))]
	end

	elem_types = split(readline(file))
	elem_nums = [parse(Int, nn) for nn in split(readline(file))]

	str = readline(file)
	if str[1] ∈ ['S', 's'] #Slective Dynamics
		str = readline(file)
	end

	if str[1] ∈ ['D', 'd'] #Direct
		location_type = "Reduced"
	elseif str[1] ∈ ['C', 'c'] #Cartesian
		location_type = "Cartesian"
	else
		error("Wrong POSCARfile.")
	end

	natoms = sum(elem_nums)
	location = Vector{Vec3{Float64}}(undef, natoms)
	for i in eachindex(location)
		location[i] = Vec3(parse.(Float64, split(readline(file))[1:3]))
	end

	close(file)

	lattice = Lattice(scale * transpose(lattice))

	name = String[]
	for (elem_type, elem_num) in zip(elem_types, elem_nums)
		append!(name, fill(elem_type, elem_num))
	end

	return Cell(lattice, location; location_type, name, period)
end
