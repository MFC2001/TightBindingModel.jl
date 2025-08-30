

function FixrExciton(bse::AbstractBSE, RSWE::RealSpaceWannierExciton, R)

	norb_ele = numorb(bse.TB)
	NR = length(R)

	location = Matrix{CartesianCoordinates{Float64}}(undef, norb_ele, NR)
	components = Matrix{ComplexF64}(undef, norb_ele, NR)

	for Ri in 1:NR, i in 1:norb_ele
		location[i, Ri] = bse.TB.lattice * (bse.TB.orb_location[i] + R[Ri])
		components[i, Ri] = RSWE(i, R[Ri], i, R[Ri])
	end

	weight = abs2.(components)

	return location, weight
end

function FixrExciton(bse::AbstractBSE, RSWE::RealSpaceWannierExciton, R; r = 0,reps = 1e-3)



	grid = Vector{CartesianCoordinates{Float64}}(undef, prod(grid_size))
	value = Vector{ComplexF64}(undef, Ngrid)

	for i in 1:Ngrid
		
	end

	r_0 = Vector{Vec3{Float64}}





	norb_ele = numorb(bse.TB)
	NR = length(R)

	location = Matrix{CartesianCoordinates{Float64}}(undef, norb_ele, NR)
	components = Matrix{ComplexF64}(undef, norb_ele, NR)

	for Ri in 1:NR, i in 1:norb_ele
		location[i, Ri] = bse.TB.lattice * (bse.TB.orb_location[i] + R[Ri])
		components[i, Ri] = RSWE(i, R[Ri], i, R[Ri])
	end

	weight = abs2.(components)

	return location, weight
end


function FixeExciton(bse::AbstractBSE, RSWE::RealSpaceWannierExciton, Rₕ::AbstractVector{<:Integer}, ei::Integer, Rₑ::AbstractVector{<:Integer})

	Rₕ = LM.gridindex(Rₕ)
	Rₕ .+= Rₑ

	norb_ele = numorb(bse.TB)
	NR = length(Rₕ)

	location = Matrix{LM.CartesianCoordinates{Float64}}(undef, norb_ele, NR)
	components = Matrix{ComplexF64}(undef, norb_ele, NR)

	for Ri in 1:NR, hi in 1:norb_ele
		location[hi, Ri] = bse.TB.lattice * (bse.TB.orb_location[hi] + Rₕ[Ri])
		components[hi, Ri] = RSWE(ei, Rₑ, hi, Rₕ[Ri])
	end

	location = vec(location)
	components = vec(components)

	weight = abs2.(components)

	return location, weight
end
function FixhExciton(bse::AbstractBSE, RSWE::RealSpaceWannierExciton, Rₑ::AbstractVector{<:Integer}, hi::Integer, Rₕ::AbstractVector{<:Integer})

	Rₑ = LM.gridindex(Rₑ)
	Rₑ .+= Rₕ

	norb_ele = numorb(bse.TB)
	NR = length(Rₑ)

	location = Matrix{LM.CartesianCoordinates{Float64}}(undef, norb_ele, NR)
	components = Matrix{ComplexF64}(undef, norb_ele, NR)

	for Ri in 1:NR, ei in 1:norb_ele
		location[ei, Ri] = bse.TB.lattice * (bse.TB.orb_location[ei] + Rₑ[Ri])
		components[ei, Ri] = RSWE(ei, Rₑ[Ri], hi, Rₕ)
	end

	location = vec(location)
	components = vec(components)

	weight = abs2.(components)

	return location, weight
end
