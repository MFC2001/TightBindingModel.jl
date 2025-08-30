struct Wannier90xsf{VT <: Number}
	lattice::Lattice{Float64}
	convvec::Lattice{Float64}
	atom_name::Vector{String}
	atom_location::Vector{ReducedCoordinates{Float64}}
	grid_size::Vec3{Int}
	origin::Vec3{Float64}
	grid_vecs::Lattice{Float64}
	value::Array{VT, 3}
end
function Readxsf(file::AbstractString, eps2::Real = 0)
	file = open(file, "r")

	#comments
	for _ in 1:4
		readline(file)
	end

	#CRYSTAL
	readline(file)

	#PRIMVEC
	readline(file)
	lattice = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		lattice[:, i] = [parse(Float64, ss) for ss in split(readline(file))]
	end
	lattice = Lattice(lattice)

	#CONVVEC
	readline(file)
	convvec = Matrix{Float64}(undef, 3, 3)
	for i in 1:3
		convvec[:, i] = [parse(Float64, ss) for ss in split(readline(file))]
	end
	convvec = Lattice(convvec)

	#PRIMCOORD
	readline(file)
	(natom, _) = parse.(Int, split(readline(file)))
	atom_name = Vector{String}(undef, natom)
	atom_location = Vector{ReducedCoordinates{Float64}}(undef, natom)
	for i in 1:natom
		line = split(readline(file))
		atom_name[i] = line[1]
		location_car = CartesianCoordinates(parse.(Float64, line[2:4]))
		atom_location[i] = lattice \ location_car
	end

	data = readdlm(file)

	close(file)

	grid_size = Int.(data[4, 1:3])
	origin = Float64.(data[5, 1:3])
	grid_vecs = Lattice(Float64.(transpose(data[6:8, 1:3])))

	value = reshape(transpose(Float64.(data[9:end-2, :])), grid_size[1], grid_size[2], grid_size[3])

	# frame = SMatrix{3, 3}(frame)
	# dvec = SMatrix{3, 3}([frame[:, 1] / Nmesh[1] frame[:, 2] / Nmesh[2] frame[:, 3] / Nmesh[3]])

	# mesh = [repeat(0:Nmesh[1]-1, outer = Nmesh[2]) repeat(0:Nmesh[2]-1, inner = Nmesh[1])]
	# mesh = [repeat(mesh, outer = (Nmesh[3], 1)) repeat(0:Nmesh[3]-1, inner = Nmesh[1] * Nmesh[2])]

	# startpoint = round.(Int, inv(dvec) * startpoint)
	# mesh = transpose(mesh) .+ startpoint

	# dV = abs(dot(dvec[:, 3], cross(dvec[:, 1], dvec[:, 2])))
	# normconst = sum(abs2, value) * dV

	# I = abs2.(value) .> eps2
	# value = value[I]
	# mesh = mesh[:, I]

	return Wannier90xsf{Float64}(lattice, convvec, atom_name, atom_location,
		grid_size, origin, grid_vecs, value)
end
