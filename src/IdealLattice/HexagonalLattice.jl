export HexagonalTB, HexagonalHR

function HexagonalTB(; a = 1, c = 20, t₀ = 0, t₁ = 1, t₂ = 0, tₕ = 0)

	l = √3 * a
	𝐚 = Vec3(l, 0, 0)
	𝐛 = l * Vec3(-1 / 2, √3 / 2, 0)
	𝐜 = Vec3(0, 0, c)

	lattice = Lattice([𝐚 𝐛 𝐜])

	location1 = Vec3(0, a, c / 2)
	location2 = [l / 2, a / 2, c / 2]

	atom_location = ReducedCoordinates.([lattice \ location1, lattice \ location2])
	orb_location = deepcopy(atom_location)

	hr = HexagonalHR(; t₀, t₁, t₂, tₕ)
	return TightBindModel(
		lattice,
		fill("C", 6),
		atom_location,
		fill("C", 6),
		orb_location,
		hr.hop,
		hr.Nhop,
		Vec3(["p", "p", "np"]))
end

function HexagonalHR(; t₀ = 0, t₁ = 1, t₂ = 0, tₕ = 0)
	path_t₀ = [
		0 0 0 1 1;
		0 0 0 2 2
	]
	value_t₀ = fill(t₀, size(path_t₀, 1))
	path_t₁ = [
		0 0 0 1 2;
		-1 0 0 1 2;
		0 1 0 1 2;
		0 0 0 2 1;
		1 0 0 2 1;
		0 -1 0 2 1
	]
	value_t₁ = fill(t₁, size(path_t₁, 1))
	path_t₂ = [
		1 1 0 1 1;
		1 0 0 1 1;
		0 1 0 1 1;
		-1 0 0 1 1;
		0 -1 0 1 1;
		-1 -1 0 1 1;
		1 1 0 2 2;
		1 0 0 2 2;
		0 1 0 2 2;
		-1 0 0 2 2;
		0 -1 0 2 2;
		-1 -1 0 2 2
	]
	value_t₂ = fill(t₂, size(path_t₂, 1))
	path_tₕ = [
		1 1 0 1 1;
		1 0 0 1 1;
		0 1 0 1 1;
		-1 0 0 1 1;
		0 -1 0 1 1;
		-1 -1 0 1 1;
		1 1 0 2 2;
		1 0 0 2 2;
		0 1 0 2 2;
		-1 0 0 2 2;
		0 -1 0 2 2;
		-1 -1 0 2 2
	]
	value_tₕ = tₕ * im * [1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1]

	path = [path_t₀; path_t₁; path_t₂; path_tₕ]
	value = [value_t₀; value_t₁; value_t₂; value_tₕ]

	I = abs.(value) .> 1e-8
	path = path[I, :]
	value = value[I]

	return HR(path, value; hrsort = "Y", buildhop = "Y")
end
