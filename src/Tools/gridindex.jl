export gridindex
"""
	gridindex(grid::AbstractVector{<:Integer}; n::Integer = 1)

	Return a grid, its datatype is like [[1,2,3],[1//2,0,0],...].
"""
function gridindex(grid::AbstractVector{<:Integer})

	start = -floor.(Int, (grid .- 1) .// 2)
	stop = ceil.(Int, (grid .- 1) .// 2)

	grid_index = [Vec3([i, j, k]) for i ∈ start[1]:stop[1], j ∈ start[2]:stop[2], k ∈ start[3]:stop[3]]

	return vec(grid_index)
end


function gridindex(cell::Cell, maxlattdist::Real)
	#Calculate the max distance, and get the logical positive lattice vectors with distances shorter than maxlattdist.
	lattice = cell.lattice
	periodicity = cell.period

	p = count(periodicity .== "p")
	if p == 3
		ucpath = gridindex(lattice, maxlattdist, Val(3))
	elseif p == 2
		if periodicity[1] == "np"
			T = 2:3
		elseif periodicity[2] == "np"
			T = [3, 1]
		elseif periodicity[3] == "np"
			T = 1:2
		else
			error("Wrong periodicity.")
		end

		ucpath = gridindex(lattice, maxlattdist, T, Val(2))

	elseif p == 1
		if periodicity[1] == "p"
			T = 1
		elseif periodicity[2] == "p"
			T = 2
		elseif periodicity[3] == "p"
			T = 3
		else
			error("Wrong periodicity.")
		end

		ucpath = gridindex(lattice, maxlattdist, T, Val(1))

	else
		error("Wrong periodicity.")
	end

	return ucpath
end
function gridindex(lattice::Lattice, maxlattdist::Real, ::Val{3})
	a₁ = lattice[:, 1]
	a₂ = lattice[:, 2]
	a₃ = lattice[:, 3]

	θ₁ = abs(calc_vecangle(a₁, a₂))
	θ₂ = abs(calc_vecangle(a₂, a₃))
	θ₃ = abs(calc_vecangle(a₃, a₁))

	dV = abs((a₁ × a₂) ⋅ a₃) / (norm(a₁) * norm(a₂) * norm(a₃))
	#It's related to the radius of biggest sphere inside parallelepiped.
	d = (2 * maxlattdist + 0.5) * maximum(sin.([θ₁, θ₂, θ₃])) / dV

	#The lattice points in a parallelepiped.
	n₁ = Int(cld(d / 2, sqrt(sum(abs2, a₁))))
	n₂ = Int(cld(d / 2, sqrt(sum(abs2, a₂))))
	n₃ = Int(cld(d / 2, sqrt(sum(abs2, a₃))))
	I = [repeat(-n₁:n₁, inner = 2 * n₂ + 1) repeat(-n₂:n₂, outer = 2 * n₁ + 1)]
	I = [repeat(I, outer = (2 * n₃ + 1, 1)) repeat(-n₃:n₃, inner = (2 * n₁ + 1) * (2 * n₂ + 1))]
	I = [SVector{3}(I[i, :]) for i in axes(I, 1)]

	D = map(x -> norm2(lattice * x), I)
	D = D .<= maxlattdist^2
	ucpath = I[D]

	return ucpath
end
function gridindex(lattice::Lattice, maxlattdist::Real, T, ::Val{2})
	lattvec1 = MMatrix{3, 3}(zeros(3, 3))
	lattvec1[:, T] .= lattice[:, T]

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
	D = D .<= maxlattdist^2
	ucpath = I[D]

	return ucpath
end
function gridindex(lattice::Lattice, maxlattdist::Real, T, ::Val{1})
	lattvec1 = MMatrix{3, 3}(zeros(3, 3))
	lattvec1[:, T] .= lattice[:, T]

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
	ucpath = I[D.<=maxlattdist^2]

	return ucpath
end
function calc_vecangle(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real})::Real
	return atan(norm(v1 × v2), v1 ⋅ v2)
end
function norm2(v::AbstractVector{<:Number})::Real
	return v ⋅ v
end
