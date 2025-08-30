function Base.show(io::IO, ::MIME"text/plain", lattice::Lattice)
	summary(io, lattice)
	println(io)
	join(io, ' ' * join(row, "  ") * '\n' for row in eachrow(parent(lattice)))
	return nothing
end
function Base.show(io::IO, ::MIME"text/plain", lattice::ReciprocalLattice)
	summary(io, lattice)
	println(io)
	join(io, ' ' * join(row, "  ") * '\n' for row in eachrow(parent(lattice)))
	return nothing
end
# function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
function Base.show(io::IO, cell::Cell)
	summary(io, cell)
	println(io)
	println(io, " lattice:")
	for row in eachrow(parent(Lattice(cell)))
		println(io, "   ", join(row, "  "))
	end
	num_atom = numatom(cell)
	println(io, " $num_atom atoms:")
	if isempty(cell.name)
		for x in cell.location
			println(io, "   ", join(x, "  "))
		end
	else
		for (name, location) in zip(cell.name, cell.location)
			println(io, "  ", name, "  ", join(location, "  "))
		end
	end
	return nothing
end