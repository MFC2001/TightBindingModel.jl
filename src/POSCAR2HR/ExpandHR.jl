
"""
	ExpandHR(unithr::HR, unitorbital::ORBITAL, unitcell::Cell, cell::Cell; keywords)

Used when cell is larger than unit cell.

Requirements:

1)`unitorbital.location` is corresponding to the orbital index in unithr one by one, of course order matters;

2)`unitorbital.belonging` is corresponding to the `unitcell`;

3)`unitcell.name` is corresponding to `unitcell.location` one by one, and this function won't use its `index`, so is `cell`;

4)`cell.name` is consistant with `unitcell.name`;

5)`lattice` of `unitcell` and `cell` need to be defined to represent the shape of cell, and all atoms should be inside the parallelepiped;

6)need to set `period` of `unitcell` and `cell` to conform to reality, 
and they should satisfy `cell.period` is equal to or less than `unitcell.period`.

Recommend that `cell` don't have lattice distortion, but `cell` can lack a few atoms.
If `cell` lack many atoms, recommend to use finduc=:custom, 
that means the unitcell is the correct unit, and don't need to rotate or translate it.


Keywords: 

`finduc` = "auto"(default), "translate" or "custom", "custom" means unitcell is the correct unitcell and don't need to rotate or translate it;

`supercell_path` = nothing(default) or a Vector, which contains some path used to expand cell for finding unitcell in cell.

`ucorientation` = (Real, Real) means angular coordinates in the spherical coordinate system, default is (0, 0), 
this function Will search the unitcell whose rotation angle is closest to this;（unusable）!

`outhr` and `outorb` = "value"(default) or "index", determine the output hr and orbital;

`outucorientation` = "N"(default) or "Y" determine if output the ucorientation or not.（unusable）!

"""
function ExpandHR(
	unithr::HR,
	unitorbital::ORBITAL,
	unitcell::Cell,
	cell::Cell; #larger than unitcell.
	finduc::AbstractString = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
	ucazimuth::Tuple{<:Real, <:Real} = (0, 0),
	outhr::AbstractString = "value",
	outorb::AbstractString = "value",
	outucorientation::AbstractString = "N",
)

	#Check the input parameters.
	outparameter = [outhr, outorb, outucorientation]
	checkKeywords(finduc, outparameter)
	(unitcell, cell) = checkCell(unitcell, cell)


	EpsPara = Dict(
		#Used by findunitcell!, permissible error between atom locations in operated unitcell and cell's supercell.
		"atomeps" => 0.1,
		#Used by findunitcell!, sum square of permissible error between atom locations in operated unitcell and cell's supercell.
		"sumeps" => 0.2,
		#Used by classifyCell, permissible error between atom locations in unitcell's supercell and cell.
		"sc_uc atom" => 0.2,
		#Used by atompath_sum in expandhr, permissible error when search cell's atom transition path in unithr.path.
		"sc_uc path" => 0.2,
		#Used by maxhrpath in hrsplit, the added value of the farthest transition lattice distance in unithr.path.
		"unithr path" => 0.1,
	)

	original_unitcell = deepcopy(unitcell)

	#This function will use name to judge whether atom is the same, so need cell.name is corresponding to unitcell.name, recommend using atom name.
	if finduc == "custom"
	else
		(unitcell, ucazimuth) = findunitcell(cell, unitcell; aimazimuth = ucazimuth, mode = finduc, supercell_path, atomeps = EpsPara["atomeps"], sumeps = EpsPara["sumeps"])
	end

	#Preprocessing unitcell HR.
	namelist!(unitcell)
	#The supercell may not be the simple repetition of unitcell, so split hr to atom pair's hrs.
	unit_atomhr = hrsplit(unithr, unitcell, unitorbital.belonging, EpsPara["unithr path"])


	#Create cell's hr from unitcell's hr.
	(hr_path, orbital_index, atom_index) = expandhr(cell, unitcell, unit_atomhr, EpsPara)


	#Generate output data.
	hr = postprocessHR(outhr, hr_path, unithr, collect(1:size(orbital_index, 1)))
	#orbital_index:[atom_index atom_orb_index unit_atom_index unit_orb_index]
	orbital = postprocessORBITAL(outorb, orbital_index, unitorbital, original_unitcell.lattice, unitcell.lattice, cell)


	if outucorientation[1] ∈ ['N', 'n']
		return hr, orbital
	elseif outucorientation[1] ∈ ['Y', 'y']
		return hr, orbital, ucorientation
	end
end

function expandhr(cell::Cell, unitcell::Cell, unit_atomhr::AtomHR, EpsPara::AbstractDict)

	#Classify atoms in cell, cell.name have been changed to namenumber.
	atom_index = classifyCell(cell, unitcell, EpsPara["sc_uc atom"])

	num_atom = length(unitcell.location)
	num_each_atom = length.(atom_index)


	#Expand supercell several times in different cases.
	maxd = maxdistance(unitcell.location)
	times = expandtimes(cell.lattice, unit_atomhr.maxlattdist + maxd, cell.period)
	#Make each atom_kind's supercell of cell.
	supercell = Vector{SuperCell}(undef, num_atom)
	for i in 1:num_atom
		#Reindex atom: 1:length[index].
		supercell[i] = SuperCell(cell.location[atom_index[i]], cell.lattice, cell.period; index = 1:num_each_atom[i])
		ExpandSuperCell!(supercell[i], times)
	end


	#Calculate some data for reindex.
	num_eachatom_orb = num_each_atom .* unit_atomhr.num_orb
	re_orb = [0; num_eachatom_orb]
	re_orb = [sum(re_orb[1:i]) for i in 1:num_atom]

	#Create cell's hr.
	hr_path = Matrix{Int}(undef, 6, 0)


	lk = ReentrantLock()
	Threads.@threads for (i, j) in [(i, j) for i in Base.OneTo(num_atom), j in Base.OneTo(num_atom)]
		if unit_atomhr.Natompath[i, j] == 0
			continue
		else
			(index1, index2) = unit_atomhr.parindex[i, j]
			dR = unitcell.location[index2] - unitcell.location[index1]
			startlocation = map(supercell_index -> cell.location[supercell_index] + dR, atom_index[i])

			pathsum = atompath_sum(startlocation, supercell[j], unit_atomhr.atompath[i, j], unitcell.lattice, unit_atomhr.maxlattdist, EpsPara["sc_uc path"])

			#Processing all orbital index.
			pathsum[4:5, :] = ((pathsum[4:5, :] .- 1) .* [unit_atomhr.num_orb[i]; unit_atomhr.num_orb[j]] .+ [re_orb[i]; re_orb[j]]) + unit_atomhr.hrorbindex[:, pathsum[6, :]]

			lock(lk) do
				hr_path = [hr_path pathsum]
			end
		end
	end


	#orbital_index:[atom_index atom_orb_index unit_atom_index unit_orb_index]
	orbital_index = Matrix{Int}(undef, 0, 4)
	for i in 1:num_atom
		unit_orb_index = unit_atomhr.orbindex[i]
		n = length(unit_orb_index)
		orbital_index = [orbital_index; repeat(atom_index[i], inner = n) repeat([collect(1:n) fill(i, n) unit_orb_index], outer = (num_each_atom[i], 1))]
	end

	return transpose(hr_path), orbital_index, atom_index
end

function maxdistance(location)
	#Return the max distance between location.
	n = length(location)

	d2 = Vector{Float64}(undef, 0)
	for i in 1:n-1
		append!(d2, map(x -> norm2(x - location[i]), location[i+1:end]))
	end

	return sqrt(maximum(d2))
end

function expandtimes(lattice, maxdist, period)

	a₁ = lattice[:, 1]
	a₂ = lattice[:, 2]
	a₃ = lattice[:, 3]

	V = abs((a₁ × a₂) ⋅ a₃)

	h₁ = V / norm(a₂ × a₃)
	h₂ = V / norm(a₃ × a₁)
	h₃ = V / norm(a₁ × a₂)

	times = Int.(cld.(maxdist, [h₁, h₂, h₃]))

	p = count(period)
	if p == 3
	elseif p == 2
		if !period[1]
			times[1] = 0
		elseif !period[2]
			times[2] = 0
		elseif !period[3]
			times[3] = 0
		end
	elseif p == 1
		if period[1]
			times[2] = 0
			times[3] = 0
		elseif period[2]
			times[1] = 0
			times[3] = 0
		elseif period[3]
			times[1] = 0
			times[2] = 0
		end
	elseif p == 0
		times .= 0
	end

	return times
end
