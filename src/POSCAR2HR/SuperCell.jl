struct SuperCell
	#Need to update.
	location::Vector{SVector{3, Float64}} #all locations in cartesian coordinate.
	path::Vector{SVector{3, Int}} #unitcell's index in supercell, corresponding to every atom
	index::Vector{Int} #supercell's atom index corresponding to ucatomindex.
	name::Vector{String} #supercell's atom name corresponding to ucatomname.
	#About the whole unitcell.
	expandtime::Base.RefValue{Int} #The times of expand = 0,1,2...
	ucpath::Vector{SVector{3, Int}} #unitcell's path in supercell

	#Don't need to update.
	lattice::Lattice{Float64}  #unitcell's lattice vector
	ucatomlocation::Vector{SVector{3, Float64}} #atom index in an unitcell
	ucatomindex::Vector{Int} #atom index in an unitcell
	ucatomname::Vector{String} #atom name in an unitcell
	period::SVector{3, Bool} #Decide how to expand supercell.
end

function SuperCell(cell::Cell)::SuperCell

	N = length(cell.location)

	return SuperCell(
		deepcopy(cell.location),
		fill(SVector{3}(0, 0, 0), N),
		deepcopy(cell.index),
		deepcopy(cell.name),
		Ref(0),
		[SVector{3}(0, 0, 0)],
		deepcopy(cell.lattice),
		deepcopy(cell.location),
		deepcopy(cell.index),
		deepcopy(cell.name),
		cell.period,
	)
end

function SuperCell(location::AbstractVector{<:AbstractVector{<:Real}},
	lattice::Lattice{<:Real},
	period::AbstractVector{<:AbstractString};
	index::AbstractVector{<:Integer} = Int[],
	name::AbstractVector{<:AbstractString} = String[],
)::SuperCell

	N = length(location)

	return SuperCell(
		deepcopy(location),
		fill(SVector{3}(0, 0, 0), N),
		deepcopy(index),
		deepcopy(name),
		Ref(0),
		[SVector{3}(0, 0, 0)],
		deepcopy(lattice),
		deepcopy(location),
		deepcopy(index),
		deepcopy(name),
		period,
	)
end


function ExpandSuperCell!(supercell::SuperCell)
	#=
	Make sure you have only used this method to expand the supercell.
	If it's periodic in three directions: expand supercell by shell, such as 1*1*1 to 3*3*3 or 3*3*3 to 5*5*5 etc.
	If it's periodic in two directions: expand supercell by circle, such as 1*1 to 3*3 or 3*3 to 5*5 etc.
	If it's periodic only in one direction: expand supercell such as 0 to -1:1 or -1:1 to -2:2.
	=#

	period = supercell.period
	p = count(period)
	if p == 3
		shell = supercell.expandtime[] + 1
		ucpath = zeros(Int, 24 * shell^2 + 2, 3)

		T = [repeat(-shell:shell, inner = 2 * shell + 1) repeat(-shell:shell, outer = 2 * shell + 1)]
		t = (2 * shell + 1)^2
		ucpath[1:t, :] = [T fill(-shell, t)]
		ucpath[t+1:2*t, :] = [T fill(shell, t)]
		n = 2 * t

		ss = shell - 1
		T = [repeat(-shell:shell, inner = 2 * shell - 1) repeat(-ss:ss, outer = 2 * shell + 1)]
		t = (2 * shell + 1) * (2 * shell - 1)
		ucpath[n+1:n+t, :] = [T[:, 1] fill(-shell, t) T[:, 2]]
		ucpath[n+t+1:n+2*t, :] = [T[:, 1] fill(shell, t) T[:, 2]]
		n = n + 2 * t

		T = [repeat(-ss:ss, inner = 2 * shell - 1) repeat(-ss:ss, outer = 2 * shell - 1)]
		t = (2 * shell - 1)^2
		ucpath[n+1:n+t, :] = [fill(-shell, t) T]
		ucpath[n+t+1:n+2*t, :] = [fill(shell, t) T]

	elseif p == 2
		circle = supercell.expandtime[] + 1
		ucpath = zeros(Int, 8 * circle, 3)
		if !period[1]
			T = 2:3
		elseif !period[2]
			T = [3, 1]
		elseif !period[3]
			T = 1:2
		else
			error("Wrong periodicity from ExpandSuperCell!.")
		end
		ucpath[1:2*circle+1, T] = [fill(-circle, 2 * circle + 1) -circle:circle]
		ucpath[2*circle+2:4*circle+2, T] = [fill(circle, 2 * circle + 1) -circle:circle]
		cc = circle - 1
		ucpath[4*circle+3:6*circle+1, T] = [-cc:cc fill(-circle, 2 * circle - 1)]
		ucpath[6*circle+2:8*circle, T] = [-cc:cc fill(circle, 2 * circle - 1)]
	elseif p == 1
		t = supercell.expandtime[] + 1
		if period[1]
			ucpath = [-t 0 0; t 0 0]
		elseif period[2]
			ucpath = [0 -t 0; 0 t 0]
		elseif period[3]
			ucpath = [0 0 -t; 0 0 t]
		else
			error("Wrong periodicity from ExpandSuperCell!.")
		end
	elseif p == 0
		return nothing
	else
		error("Wrong periodicity from expandSuperCell!")
	end

	ucpath = [SVector{3}(ucpath[i, :]) for i in axes(ucpath, 1)]
	ExpandSuperCell!(supercell, ucpath)

	return nothing
end

function ExpandSuperCell!(supercell::SuperCell, times::AbstractVector{<:Integer})

	if any(times .< 0)
		error("Please use right parameters!")
	end

	aimpath = [repeat(-times[1]:times[1], inner = 2 * times[2] + 1) repeat(-times[2]:times[2], outer = 2 * times[1] + 1)]
	aimpath = [repeat(aimpath, inner = (2 * times[3] + 1, 1)) repeat(-times[3]:times[3], outer = size(aimpath, 1))]
	aimpath = [SVector{3}(aimpath[i, :]) for i in axes(aimpath, 1)]


	ExpandSuperCell!(supercell, aimpath)

	return nothing
end

function ExpandSuperCell!(supercell::SuperCell, aimlocation::AbstractVector{<:Real})
	#=
	Will try to expand supercell to include aim atom.
	=#

	dR = map(x -> supercell.lattice \ (aimlocation - x), supercell.ucatomlocation)
	aimI = map(x -> round.(x), dR)

	(_, i) = findmin(i -> norm(aimI[i] - dR[i]), eachindex(dR))
	aimI = Int.(aimI[i])

	period = supercell.period
	p = count(period)
	if p == 3
		T = [repeat(-1:1, inner = 3) repeat(-1:1, outer = 3)]
		T = [repeat(T, inner = (3, 1)) repeat(-1:1, outer = 9)]

	elseif p == 2
		if !period[1]
			T = [zeros(Int, 9) repeat(-1:1, inner = 3) repeat(-1:1, outer = 3)]
		elseif !period[2]
			T = [repeat(-1:1, inner = 3) zeros(Int, 9) repeat(-1:1, outer = 3)]
		elseif !period[3]
			T = [repeat(-1:1, inner = 3) repeat(-1:1, outer = 3) zeros(Int, 9)]
		else
			error("Wrong periodicity from ExpandSuperCell!.")
		end
	elseif p == 1
		if period[1]
			T = [-1:1 zeros(Int, 3, 2)]
		elseif period[2]
			T = [zeros(Int, 3) -1:1 zeros(Int, 3)]
		elseif period[3]
			T = [zeros(Int, 3, 2) -1:1]
		else
			error("Wrong periodicity from ExpandSuperCell!.")
		end
	elseif p == 0
		return nothing
	else
		error("Wrong periodicity from expandSuperCell!")
	end
	aimpath = T .+ transpose(aimI)
	aimpath = [SVector{3}(aimpath[i, :]) for i in axes(aimpath, 1)]


	if isempty(aimpath)
		return nothing
	else
		ExpandSuperCell!(supercell, aimpath)
		return nothing
	end
end

function ExpandSuperCell!(supercell::SuperCell, ucpath::AbstractVector{<:SVector{3, <:Integer}})
	#ucpath = supercell's unitcell index that is needed to add.

	ucpath = setdiff(ucpath, supercell.ucpath)

	ucatomnum = length(supercell.ucatomlocation)
	I = Base.OneTo(ucatomnum)

	for i in eachindex(ucpath)
		L = supercell.lattice * ucpath[i]
		append!(supercell.location, supercell.location[I] .+ [L])
		append!(supercell.path, fill(ucpath[i], ucatomnum))
		append!(supercell.index, supercell.ucatomindex)
		append!(supercell.name, supercell.ucatomname)
	end

	supercell.expandtime[] += 1
	append!(supercell.ucpath, ucpath)

	return nothing
end
