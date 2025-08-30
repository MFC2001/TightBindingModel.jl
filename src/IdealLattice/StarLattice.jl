export StarTB, StarCell, StarHR

function StarTB(; a = 1, b = 2, c = 20,
	t₀ = 0, t₁ = 1.5, t₂ = 1, t₃ = 0, t₄ = 0, t₅ = 0, t_inverse = 0)

	cell = StarCell(; a, b, c)
	orb_location = deepcopy(cell.location)
	orb_name = deepcopy(cell.name)

	hr = StarHR(; t₀, t₁, t₂, t₃, t₄, t₅, t_inverse)

	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orb_name,
		orb_location,
		HermitianReciprocalHoppings(hr),
		cell.period)
end
function StarCell(; a = 1, b = 2, c = 20)

	l = √3 * a + 2 * b
	𝐚 = Vec3(l, 0, 0)
	𝐛 = l * Vec3(-1 / 2, √3 / 2, 0)
	𝐜 = Vec3(0, 0, c)

	lattice = Lattice([𝐚 𝐛 𝐜])

	location1 = Vec3(l / 2, a / 2, c / 2)
	location2 = location1 + [b / 2, √3 * b / 2, 0]
	location3 = location1 + [-b / 2, √3 * b / 2, 0]

	location5 = 𝐛 / 2 + a / 2 * [√3 / 2, 1 / 2, 0] + 𝐜 / 2
	location4 = location5 + [b, 0, 0]
	location6 = location4 + [-b / 2, √3 * b / 2, 0]

	location = ReducedCoordinates.([lattice \ location1, lattice \ location2, lattice \ location3, lattice \ location4, lattice \ location5, lattice \ location6])

	return Cell(lattice, location; location_type = "ReducedCoordinates", name = fill("X", 6), period = Bool[1, 1, 0])
end
function StarHR(; t₀ = 0, t₁ = 1.5, t₂ = 1, t₃ = 0, t₄ = 0, t₅ = 0, t_inverse = 0)
	path_t₀ = [
		0 0 0 1 1;
		0 0 0 2 2;
		0 0 0 3 3;
		0 0 0 4 4;
		0 0 0 5 5;
		0 0 0 6 6
	]
	value_t₀ = fill(t₀, size(path_t₀, 1))
	path_t₁ = [
		0 0 0 3 4;
		0 0 0 4 3;
		0 -1 0 1 6;
		1 0 0 2 5;
		-1 0 0 5 2;
		0 1 0 6 1
	]
	value_t₁ = fill(t₁, size(path_t₁, 1))
	path_t₂ = [
		0 0 0 1 2;
		0 0 0 2 1;
		0 0 0 1 3;
		0 0 0 3 1;
		0 0 0 2 3;
		0 0 0 3 2;
		0 0 0 4 5;
		0 0 0 5 4;
		0 0 0 4 6;
		0 0 0 6 4;
		0 0 0 5 6;
		0 0 0 6 5
	]
	value_t₂ = fill(t₂, size(path_t₂, 1))
	path_t₃ = [
		0 0 0 1 4;
		1 0 0 1 5;
		0 -1 0 1 4;
		0 -1 0 1 5;
		0 0 0 2 4;
		0 -1 0 2 6;
		1 0 0 2 4;
		1 0 0 2 6;
		0 0 0 3 5;
		0 0 0 3 6;
		1 0 0 3 5;
		0 -1 0 3 6;
		0 0 0 4 1;
		0 0 0 4 2;
		0 1 0 4 1;
		-1 0 0 4 2;
		0 0 0 5 3;
		0 1 0 5 1;
		-1 0 0 5 1;
		-1 0 0 5 3;
		0 0 0 6 3;
		-1 0 0 6 2;
		0 1 0 6 2;
		0 1 0 6 3
	]
	value_t₃ = fill(t₃, size(path_t₃, 1))
	path_t₄ = [
		0 0 0 1 5;
		1 0 0 1 4;
		0 0 0 2 6;
		0 -1 0 2 4;
		0 -1 0 3 5;
		1 0 0 3 6;
		-1 0 0 4 1;
		0 1 0 4 2;
		0 0 0 5 1;
		0 1 0 5 3;
		0 0 0 6 2;
		-1 0 0 6 3
	]
	value_t₄ = fill(t₄, size(path_t₄, 1))
	path_t₅ = [
		-1 -1 0 1 2;
		0 -1 0 1 3;
		1 1 0 2 1;
		1 0 0 2 3;
		0 1 0 3 1;
		-1 0 0 3 2;
		1 0 0 4 5;
		0 -1 0 4 6;
		-1 0 0 5 4;
		-1 -1 0 5 6;
		0 1 0 6 4;
		1 1 0 6 5
	]
	value_t₅ = fill(t₅, size(path_t₅, 1))

	path_t_inverse = [
		0 0 0 1 1;
		0 0 0 2 2;
		0 0 0 3 3;
		0 0 0 4 4;
		0 0 0 5 5;
		0 0 0 6 6
	]
	value_t_inverse = [t_inverse, t_inverse, t_inverse, -t_inverse, -t_inverse, -t_inverse]

	path = [path_t₀; path_t₁; path_t₂; path_t₃; path_t₄; path_t₅; path_t_inverse]
	value = [value_t₀; value_t₁; value_t₂; value_t₃; value_t₄; value_t₅; value_t_inverse]

	I = abs.(value) .> 1e-8
	path = path[I, :]
	value = value[I]

	return HR(path, value; orbindex = collect(1:6), hrsort = "Y")
end
function StarHR(::Val{:Haldane}, t_uu::Real,t_dd::Real = -t_uu)
	haldane_uu_path = [
		0 0 0 1 4 1;
		1 0 0 1 5 -1;
		0 -1 0 1 4 1;
		0 -1 0 1 5 -1;
		0 0 0 2 4 -1;
		0 -1 0 2 6 1;
		1 0 0 2 4 -1;
		1 0 0 2 6 1;
		0 0 0 3 5 1;
		0 0 0 3 6 -1;
		1 0 0 3 5 1;
		0 -1 0 3 6 -1;
		0 0 0 4 1 -1;
		0 0 0 4 2 1;
		0 1 0 4 1 -1;
		-1 0 0 4 2 1;
		0 1 0 5 1 1;
		0 0 0 5 3 -1;
		-1 0 0 5 1 1;
		-1 0 0 5 3 -1;
		0 -1 0 6 2 -1;
		0 0 0 6 3 1;
		0 1 0 6 2 -1;
		0 1 0 6 3 1
	]
	haldane_uu_value = haldane_uu_path[:, 6] .* (complex(0, t_uu))
	haldane_dd_value = haldane_uu_path[:, 6] .* (complex(0, t_dd))
	haldane_uu_path = haldane_uu_path[:, 1:5]
	haldane_dd_path = copy(haldane_uu_path)
	haldane_dd_path[:, 4:5] .+= 6

	path = [haldane_uu_path; haldane_dd_path]
	value = [haldane_uu_value; haldane_dd_value]

	return HR(path, value;orbindex = collect(1:12))
end
