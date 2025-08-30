"""
	findunitcell(cell::Cell, unitcell::Cell; aimazimuth::Tuple{<:Real, <:Real} = (0, 0), mode = "auto", atomeps::Real = 0.1, sumeps::Real = 0)

	In any supercell, find an unitcell. 
	Will search the unitcell whose rotation angle is closest to aimazimuth.(unusable)

	`mode` = "auto", "translate"

"""
function findunitcell(cell::Cell, unitcell::Cell;
	aimazimuth::Tuple{<:Real, <:Real} = (0, 0),
	mode = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
	atomeps::Real = 0.1,
	sumeps::Real = 0)

	n = length(unitcell.location)

	#TODO These parameters is to be continued.
	if n < 4
		t = atomeps
	else
		t = n * atomeps / 5
	end
	if sumeps < t
		sumeps = t
	end

	#unitcell should have 4 atoms at least.
	if n == 1
		expandcell = expanduc(unitcell, 3)
	elseif n < 4
		expandcell = expanduc(unitcell, 1)
	else
		expandcell = unitcell
	end

	supercell = SuperCell(cell)
	if isnothing(supercell_path)
		ExpandSuperCell!(supercell)
	else
		try
			supercell_path = map(SVector{3, Int}, supercell_path)
		catch e
			error("Make sure supercell_path is a vector of vector with 3 integers!")
		end
		ExpandSuperCell!(supercell, supercell_path)
	end


	mode = Symbol("finduc", mode)
	(center, azimuth) = eval(mode)(supercell, expandcell, aimazimuth, atomeps, sumeps)

	rm = azimuth[2]
	location = map(x -> CartesianCoordinates(rm * (x - unitcell.location[1]) + center), unitcell.location)
	lattice = Lattice(rm * parent(unitcell.lattice))
	unitcell = Cell(lattice, location, unitcell.name, unitcell.index, unitcell.period)

	return unitcell, azimuth[1]
end

function finducauto(supercell::SuperCell, expandcell::Cell, aimazimuth::Tuple{Real, Real}, atomeps::Real, sumeps::Real)

	N = length(expandcell.location)
	atomeps = atomeps^2

	sclocat = supercell.location
	scname = supercell.name

	#Some value will be used later.
	dlocation = map(x -> x - expandcell.location[1], expandcell.location)
	d₂ = sum(abs2, dlocation[2])
	D₂ = dlocation[2] / sqrt(d₂)
	d₃ = sum(abs2, dlocation[3])
	Δ = dlocation[3] - (dlocation[3] ⋅ D₂) * D₂
	δ = sum(abs2, Δ)
	I₁ = findall(scname .== expandcell.name[1])
	I₂ = findall(scname .== expandcell.name[2])
	I₃ = findall(scname .== expandcell.name[3])

	center = nothing
	azimuth = Vector{Pair}(undef, 0)
	for centerindex in I₁
		center = sclocat[centerindex]

		#Find samely separated atom with same expandcell.name[2]
		D = map(x -> norm2(x - center), sclocat[I₂])
		II₂ = I₂[abs.(D .- d₂).<0.1]
		#Find samely separated atom with same expandcell.name[3]
		T = map(x -> x - center, sclocat[I₃])
		D = map(norm2, T)
		I = abs.(D .- d₃) .< 0.1
		T = T[I]
		II₃ = I₃[I]

		#Find valid azimuth angle (θ, ϕ) and rotate angle ω.
		for i₂ in II₂
			k = sclocat[i₂] - center
			k = k / norm(k)
			rm1 = RM(D₂, k)
			D₃ = rm1 * Δ
			T₃ = k × D₃

			#Find properly atom3.
			TT = map(x -> x - (x ⋅ k) * k, T)
			D = map(norm2, TT)
			I = abs.(D .- δ) .< 0.1
			TT = TT[I]
			III₃ = II₃[I]

			for i₃ in eachindex(III₃)
				ω = atan(TT[i₃] ⋅ T₃, TT[i₃] ⋅ D₃)
				rm2 = RM(k, ω)

				rm = rm2 * rm1
				rotat_dlocation = map(x -> rm * x + center, dlocation)

				sum_d = 0
				for j in 4:N
					I = scname .== expandcell.name[j]
					D = map(x -> norm2(x - rotat_dlocation[j]), sclocat[I])
					d = minimum(D)
					if d < atomeps
						sum_d += sqrt(d)
					else
						sum_d = sumeps + 1
						break
					end
				end

				if sum_d < sumeps
					#TODO How to represent these three angles, Euler angles? or others.
					(θ, ϕ) = coord2azimuth(rm * [0, 0, 1])
					push!(azimuth, Pair([θ, ϕ, ω], rm))
				end
			end
		end

		if !(isempty(azimuth))
			break
		end
	end

	#Find the direction closest to the specified direction.
	# azimuth = findcloset(azimuth, aimazimuth)
	if isempty(azimuth)
		error("Cannot find any unitcell in supercell, they may be incompatible.")
	end

	return center, azimuth[1]
end
function coord2azimuth(k::AbstractVector{<:Real})
	x = k[1]
	y = k[2]
	z = k[3]
	θ = atan(sqrt(x^2 + y^2), z)
	ϕ = atan(y, x)
	return θ, ϕ
end


function finductranslate(supercell::SuperCell, expandcell::Cell, aimazimuth::Tuple{Real, Real}, atomeps::Real, sumeps::Real)

	N = length(expandcell.location)
	atomeps = atomeps^2

	#Some value that will be used later.
	dlocation = map(x -> x - expandcell.location[1], expandcell.location)

	finduc = "N"
	center = nothing
	for centerindex in findall(supercell.name .== expandcell.name[1])
		center = supercell.location[centerindex]

		translate_location = map(x -> x + center, dlocation)

		sum_d = 0
		for j in 2:N
			T = supercell.name .== expandcell.name[j]
			D = map(norm2, supercell.location[T] .- [translate_location[j]])
			d = minimum(D)
			if d < atomeps
				sum_d += sqrt(d)
			else
				sum_d = sumeps + 1
				break
			end
		end

		if sum_d < sumeps
			finduc = "Y"
			break
		end
	end

	if finduc == "N"
		error("Cannot find unit cell with mode=\"translate\"!")
	end


	return center, Pair([0, 0, 0], [1 0 0; 0 1 0; 0 0 1])
end

function expanduc(cell::Cell, n::Integer)

	period = cell.period
	p = count(period)
	if p == 3
		T = [
			1 0 0;
			0 1 0;
			1 1 0
		]

	elseif p == 2
		if !period[1]
			T = [0 1 0; 0 0 1; 0 1 1]
		elseif !period[2]
			T = [1 0 0; 0 0 1; 1 0 1]
		elseif !period[3]
			T = [1 0 0; 0 1 0; 1 1 0]
		end
	elseif p == 1
		if period[1]
			T = [1 0 0; -1 0 0; 2 0 0; -2 0 0]
		elseif period[2]
			T = [0 1 0; 0 -1 0; 0 2 0; 0 -2 0]
		elseif period[3]
			T = [0 0 1; 0 0 -1; 0 0 2; 0 0 -2]
		end
	elseif p == 0
		error("Wrong periodicity of unitcell, it should be periodic in at least one dimension.")
	end

	#Make sure the first part of location is original unitcell.
	location = deepcopy(cell.location)
	for i in 1:n
		append!(location, map(x -> x + cell.lattice * T[i, :], location))
	end
	name = repeat(cell.name, n + 1)
	index = repeat(cell.index, n + 1)

	return Cell(cell.lattice, location; name, index, period)
end

"""
	RM(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::AbstractMatrix

	For x,y ∈ Rⁿ and |x|=|y|=1, y=R(x,y)x.
"""
function RM(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
	xy = x ⋅ y
	if abs(xy + 1) < 1e-12
		local a
		for _ in 1:100
			a = rand(eltype(x), length(x))
			if abs(a ⋅ x) < 1
				break
			end
		end
		n = a - (a ⋅ x) * x
		n = n / norm(n)
		return RM(n, π)
	else
		T = y * transpose(x) - x * transpose(y)
		return I + T + T^2 ./ (1 + x ⋅ y)
	end
end
"""
	RM(n::AbstractVector{<:Real}, ω::Real)::AbstractMatrix{3,3}

	e^{-iω n⋅T}
"""
function RM(n::AbstractVector{<:Real}, ω::Real)
	iT₁ = zeros(Int, 3, 3)
	iT₁[2, 3] = 1
	iT₁[3, 2] = -1
	iT₂ = zeros(Int, 3, 3)
	iT₂[1, 3] = -1
	iT₂[3, 1] = 1
	iT₃ = zeros(Int, 3, 3)
	iT₃[1, 2] = 1
	iT₃[2, 1] = -1
	return exp(-ω * (n[1] * iT₁ + n[2] * iT₂ + n[3] * iT₃))
end

