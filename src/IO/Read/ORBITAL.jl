export ReadORBITAL

"""
	ReadORBITAL(ORBITALfile::AbstractString)::ORBITAL

Read wannier90_centres.xyz, return atoms & locations corresponding to each orbital.
"""
function ReadORBITAL(file::AbstractString)::ORBITAL

	file = open(file, "r")
	N = parse(Int, readline(file))
	readline(file)

	locations = readdlm(file)
	close(file)

	I = occursin.(r"^X", locations[:, 1])
	orbital_index = locations[I, 1]
	orbital_location = convert.(Float64, locations[I, 2:4])
	atom_name = convert.(String, locations[.!I, 1])
	atom_location = convert.(Float64, locations[.!I, 2:4])

	orbital_location = [CartesianCoordinates(r) for r in eachrow(orbital_location)]
	atom_location = [CartesianCoordinates(r) for r in eachrow(atom_location)]
	

	norb = count(I)
	name = Vector{String}(undef, norb)
	belonging = Vector{Int}(undef, norb)
	for i in 1:norb
		(_, I) = findmin(x -> sum(abs2, x - orbital_location[i]), atom_location)
		name[i] = atom_name[I]
		belonging[i] = I
	end

	#Get orbital index from "X*"
	index = match.(r"^X([0-9]+)", orbital_index)
	if count(index .== nothing) > 0
		orbital_index = collect(1:norb)
	else
		for i in eachindex(index)
			orbital_index[i] = index[i].captures[1]
		end
		orbital_index = parse.(Int, orbital_index)
	end

	return ORBITAL(orbital_location; name, index = orbital_index, atom_location, atom_name, belonging)
end
