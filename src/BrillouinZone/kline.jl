export kline
"""
	kline(nk::Integer, cell)
	kline(point::AbstractVector{<:Pair}, nk::Integer, cell::Cell)

`point` = [ "Name" => [frac_coord],...], [frac_coord] is a vector with 3 elements.
"""
function kline(nk::Integer, TB::AbstractTightBindModel)
	return kline(nk, TB.lattice, TB.period)
end
function kline(nk::Integer, cell::Cell)
	return kline(nk, cell.lattice, cell.period)
end
function kline(nk::Integer, lattice::Lattice, period)

	𝐚, 𝐛, 𝐜 = basisvectors(lattice)

	p = count(period)

	if p == 0
		point = Vector{Pair}(undef, 0)
	elseif p == 1
		if period[1]
			name = "X"
			point = [
				"-" * name => [-0.5, 0, 0],
				"Γ" => [0, 0, 0],
				name => [0.5, 0, 0],
			]
		elseif period[2]
			name = "Y"
			point = [
				"-" * name => [0, -0.5, 0],
				"Γ" => [0, 0, 0],
				name => [0, 0.5, 0],
			]
		elseif period[3]
			name = "Z"
			point = [
				"-" * name => [0, 0, -0.5],
				"Γ" => [0, 0, 0],
				name => [0, 0, 0.5],
			]
		else
			error("Wrong period from kline.")
		end
	elseif p == 2
		if !period[1]
			a₁ = 𝐛
			a₂ = 𝐜
		elseif !period[2]
			a₁ = 𝐚
			a₂ = 𝐜
		elseif !period[3]
			a₁ = 𝐚
			a₂ = 𝐛
		else
			error("Wrong period from kline.")
		end

		point = Vector{Pair}(undef, 0)
		if !(isHexagon!(point, a₁, a₂) || isSquare!(point, a₁, a₂))
			point = [
				"Γ" => [0, 0],
				"X" => [0.5, 0],
				"M" => [0.5, 0.5],
				"Y" => [0, 0.5],
				"Γ" => [0, 0],
				"M" => [0.5, 0.5],
			]
		end
		point = expandpoint(point, period)
	elseif p == 3
		error("To be continued.")
	else
		error("Wrong period from kline.")
	end

	return kline(point, nk, lattice)
end
function kline(point::AbstractVector{<:Pair}, nk::Integer, TB::AbstractTightBindModel)
	return kline(point, nk, TB.lattice)
end
function kline(point::AbstractVector{<:Pair}, nk::Integer, cell::Cell)
	return kline(point, nk, cell.lattice)
end
function kline(point::AbstractVector{<:Pair}, nk::Integer, lattice::Lattice)

	reciprocallattice = reciprocal(lattice)
	basis = reciprocallattice.data

	# line is Vector{Float64}
	if length(point) == 0
		return Kline(; basis, kdirect = [Vec3(0, 0, 0)], line = [0.0], name = ["Γ"], index = [1])
	elseif length(point) == 1
		return Kline(; basis, kdirect = Vec3(point[1][2]), line = [0.0], name = [point[1][1]], index = [1])
	end

	dk = norm(basis * (point[2][2] - point[1][2])) / nk

	point_coord = similar(point, Vec3{Float64})
	point_name = similar(point, String)
	for i in eachindex(point)
		point_coord[i] = Vec3(basis * point[i][2])
		point_name[i] = point[i][1]
	end
	point_frac = [point[i][2] for i in eachindex(point)]


	kdirect, line, index = kpoint2kline(point_coord, point_frac, dk)

	return Kline(; basis, kdirect, line, name = point_name, index)
end
function expandpoint(point, periodicity)

	l = length(point[1].second)
	p = count(periodicity)

	if l ≠ p
		error("Wrong kpoint, please use the same dimension with periodicity (with the addition of surfdirection for surfkline).")
	end

	if p == 3
		return point
	elseif p == 2
		if !periodicity[1]
			T = 2:3
		elseif !periodicity[2]
			#Watch out!
			T = [1, 3]
		elseif !periodicity[3]
			T = 1:2
		else
			error("Wrong periodicity from kline.")
		end
	elseif p == 1
		if periodicity[1]
			T = [1]
		elseif periodicity[2]
			T = [2]
		elseif periodicity[3]
			T = [3]
		else
			error("Wrong periodicity from kline.")
		end
	elseif p == 0
		error("Please don't input kpoint for cluster.")
	end

	newpoint = similar(point, Pair{String, Vector{Float64}})
	for i in eachindex(point)
		ll = length(point[i].second)
		if ll ≠ l
			error("Wrong kpoint, please use the same dimension.")
		elseif ll == 3
			continue
		end

		newpointi = zeros(3)
		newpointi[T] .= point[i].second

		newpoint[i] = point[i].first => Vec3(newpointi)
	end

	return newpoint
end
function isSquare!(point, a₁, a₂)

	if abs(a₁ ⋅ a₂) < 1e-4 && norm(a₁) == norm(a₂)
		append!(point, [
			"Γ" => [0, 0],
			"X" => [0.5, 0],
			"M" => [0.5, 0.5],
			"Γ" => [0, 0],
		])
		return true
	end

	return false
end
function isHexagon!(point, a₁, a₂)

	A₁ = norm(a₁)
	A₂ = norm(a₂)
	if abs(A₁ - A₂) < 1e-5
		θ = acos((a₁ ⋅ a₂) / (A₁ * A₂))
		if abs(θ - π / 3) < 1e-3
			append!(point, [
				"Γ" => [0, 0],
				"M" => [0.5, 0.5],
				"K" => [2 / 3, 1 / 3],
				"Γ" => [0, 0],
			])
			return true
		elseif abs(θ - 2 * π / 3) < 1e-3
			append!(point, [
				"Γ" => [0, 0],
				"M" => [0.5, 0],
				"K" => [1 / 3, 1 / 3],
				"Γ" => [0, 0],
			])
			return true
		end
	end

	return false
end

function kpoint2kline(point, point_frac, dk)

	m = length(point)
	index = zeros(Int, m)

	#Sequentially generating other points.
	for i in 2:m
		d = norm(point[i] - point[i-1])
		index[i] = round(Integer, d / dk)
		if index[i] == 0
			index[i] = 1
		end
	end

	#Total number of points = Total segments + 1.
	N = sum(index) + 1
	#Making kline coordinates & index of High-symmetry points(kn).
	kx = zeros(N)
	ky = zeros(N)
	kz = zeros(N)
	kx_frac = zeros(N)
	ky_frac = zeros(N)
	kz_frac = zeros(N)
	#Starting point.
	t = 1
	index[1] = t
	#Sequentially generating other points.
	#tt is the number of segments about each line, t is previous end point.
	#Overwrite High-symmetry point repeatedly, so t:t+tt not t+1:t+tt.
	for i ∈ 2:m
		tt = index[i]
		kx[t:t+tt] = range(point[i-1][1], point[i][1], length = tt + 1)
		ky[t:t+tt] = range(point[i-1][2], point[i][2], length = tt + 1)
		kz[t:t+tt] = range(point[i-1][3], point[i][3], length = tt + 1)
		kx_frac[t:t+tt] = range(point_frac[i-1][1], point_frac[i][1], length = tt + 1)
		ky_frac[t:t+tt] = range(point_frac[i-1][2], point_frac[i][2], length = tt + 1)
		kz_frac[t:t+tt] = range(point_frac[i-1][3], point_frac[i][3], length = tt + 1)
		t = t + tt
		index[i] = t
	end
	k = [Vec3(k) for k in eachrow([kx ky kz])]

	line = zeros(N)
	for i in 2:N
		d = norm(k[i] - k[i-1])
		line[i] = line[i-1] + d
	end

	k_frac = [ReducedCoordinates(k) for k in eachrow([kx_frac ky_frac kz_frac])]

	return k_frac, line, index
end
