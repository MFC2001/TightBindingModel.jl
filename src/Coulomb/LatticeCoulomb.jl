export LatticeCoulomb
"""
	LatticeCoulomb(Coulomb::Function, kgrid::MonkhorstPack, TB::AbstractTightBindModel)::HR
	LatticeCoulomb(Coulomb::Function, kgrid::MonkhorstPack, cell::Cell,	orblocation::AbstractVector)::HR

"""
function LatticeCoulomb(Coulomb::Function, kgrid::MonkhorstPack, TB::AbstractTightBindModel)::HR

	ucpath = gridindex(kgrid.kgrid_size; n = 1)

	return LatticeCoulomb(Coulomb, TB.lattice, ucpath, TB.orb_location)
end
function LatticeCoulomb(Coulomb::Function, kgrid::MonkhorstPack, cell::Cell,
	orblocation::AbstractVector)::HR

	ucpath = gridindex(kgrid.kgrid_size; n = 1)

	return LatticeCoulomb(Coulomb, cell.lattice, ucpath, orblocation)
end
"""
	LatticeCoulomb(Coulomb::Function, lattice::Lattice, ucpath::AbstractVector, orblocation::AbstractVector; n::Integer = 1)::HR

	Coulomb(r) = V(r).

	Using the positive lattice vectors, return a interaction HR.
"""
function LatticeCoulomb(
	Coulomb::Function,
	lattice::Lattice,
	ucpath::AbstractVector{<:AbstractVector{<:Integer}},
	orblocation::AbstractVector{<:ReducedCoordinates{<:Real}},
)::HR

	orblocation = map(x -> lattice * x, orblocation)

	return LatticeCoulomb(Coulomb, lattice, ucpath, orblocation)
end
function LatticeCoulomb(
	Coulomb::Function,
	lattice::Lattice,
	ucpath::AbstractVector{<:AbstractVector{<:Integer}},
	orblocation::AbstractVector{<:CartesianCoordinates{<:Real}},
)::HR

	norb = length(orblocation)
	N = size(ucpath, 1)

	R = map(x -> lattice * x, ucpath)

	B = Matrix{SVector{3, Float64}}(undef, norb, norb)
	for i in 1:norb, j in 1:i
		dorb = orblocation[j] - orblocation[i]
		B[i, j] = dorb
		B[j, i] = -dorb
	end

	Uhrpath = Matrix{Int}(undef, 5, N * norb^2)
	Uhrvalue = Vector{Float64}(undef, N * norb^2)

	orbpath = [repeat(1:norb, inner = norb) repeat(1:norb, outer = norb)]'
	Threads.@threads for i in 1:N
		I = norb^2*(i-1)+1:norb^2*i
		Uhrpath[:, I] = [repeat(ucpath[i], inner = (1, norb^2)); orbpath]
		Uhrvalue[I] = map(x -> Coulomb(B[x[1], x[2]] + R[i]), eachcol(orbpath))
	end

	return HR(Uhrpath', Uhrvalue)
end
