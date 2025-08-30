
"""
	ContractHR(superhr::HR, superorbital::ORBITAL, supercell::Cell, cell::Cell; keywords)

Used when cell is equal to or less than supercell.

Requirements:

1)`superorbital.location` is corresponding to the orbital index in superhr one by one, of course order matters;

2)`superorbital.belonging` is corresponding to the `supercell`;

3)`supercell.name` is corresponding to `supercell.location` one by one, and this function won't use its `index`, so is `cell`;

4)`cell.name` is consistant with `supercell.name`;

5)`lattice` of `supercell` and `cell` need to be defined to represent the shape of cell, and all atoms should be inside the parallelepiped;

6)need to set `period` of `supercell` and `cell` to conform to reality, 
and they should satisfy `cell.period` is equal to or less than `supercell.period`.

Recommend that `cell` don't have lattice distortion, but `cell` can lack a few atoms.
If `cell` lack many atoms, recommend to use finduc=:custom, 
that means we don't need to rotate or translate the cell.


Keywords: 

`finduc` = "auto"(default), "translate" or "custom", "custom" means we don't need to rotate or translate the cell;

`supercell_path` = nothing(default) or a Vector, which contains some path used to expand supercell for finding cell in supercell.

`ucorientation` = (Real, Real) means angular coordinates in the spherical coordinate system, default is (0, 0), 
this function Will search the unitcell whose rotation angle is closest to this;（unusable）!

`outhr` and `outorb` = "value"(default) or "index", determine the output hr and orbital;

`outucorientation` = "N"(default) or "Y" determine if output the ucorientation or not.（unusable）!

"""
function ContractHR(
	superhr::HR,
	superorbital::ORBITAL,
	supercell::Cell,
	cell::Cell; #equal to or less than supercell.
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
	(supercell, cell) = checkCell(supercell, cell)


	EpsPara = Dict(
		#Used by findunitcell!, permissible error between atom locations in operated unitcell and cell's supercell.
		"atomeps" => 0.1,
		#Used by findunitcell!, sum square of permissible error between atom locations in operated unitcell and cell's supercell.
		"sumeps" => 0.2,
		#Used by classifyCell, permissible error between atom locations in unitcell's supercell and cell.
		"sc_uc atom" => 0.2,
		#Used by atompath_sum in poscar2hr, permissible error when search cell's atom transition path in unithr.path.
		"sc_uc path" => 0.2,
		#Used by maxhrpath in hrsplit, the added value of the farthest transition lattice distance in unithr.path.
		"unithr path" => 0.1,
	)

	# Need cell without change to generate orbital.
	original_cell = deepcopy(cell)

	#This function will use name to judge whether atom is the same, so need cell.name is corresponding to unitcell.name, recommend using atom name.
	if finduc == "custom"
	else
		(cell, ucazimuth) = findunitcell(supercell, cell; aimazimuth = ucazimuth, mode = finduc, supercell_path, atomeps = EpsPara["atomeps"], sumeps = EpsPara["sumeps"])
	end

	#Preprocessing unitcell HR.
	namelist!(supercell)
	#The supercell may not be the simple repetition of unitcell, so split hr to atom pair's hrs.
	super_atomhr = hrsplit(superhr, supercell, superorbital.belonging, EpsPara["unithr path"])


	#Create cell's hr from unitcell's hr.
	(hr_path, orbital_index, atom_index) = contracthr(cell, supercell, super_atomhr, EpsPara)


	#Generate output data.
	hr = postprocessHR(outhr, hr_path, superhr, collect(1:size(orbital_index, 1)))
	#orbital_index:[atom_index atom_orb_index super_atom_index super_orb_index]
	orbital = postprocessORBITAL(outorb, orbital_index, superorbital, cell.lattice, original_cell.lattice, original_cell)


	if outucorientation[1] ∈ ['N', 'n']
		return hr, orbital
	elseif outucorientation[1] ∈ ['Y', 'y']
		return hr, orbital, ucorientation
	end
end

function contracthr(cell::Cell, supercell::Cell, super_atomhr::AtomHR, EpsPara::AbstractDict)

	num_atom = length(cell.location)
	super_num_atom = length(supercell.location)

	#Classify atoms in cell, cell.name have been changed to namenumber.
	atom_index = classifyCell(supercell, cell, EpsPara["sc_uc atom"])

	atom_index_super = Vector{Int}(undef, super_num_atom)
	for i in eachindex(atom_index)
		for j in atom_index[i]
			atom_index_super[j] = i
		end
	end

	TM = inv(cell.lattice.data) * supercell.lattice.data

	atompath = [Vector{AtomPath}(undef, 0) for _ in 1:num_atom, _ in 1:num_atom]
	lk = ReentrantLock()
	Threads.@threads for (i, j_super) in [(i, j) for i in Base.OneTo(num_atom), j in Base.OneTo(super_num_atom)]
		i_super = atom_index[i][1]
		j = atom_index_super[j_super]

		dR_super = supercell.location[j_super] - supercell.location[i_super]
		dR = cell.location[j] - cell.location[i]
		Δ = cell.lattice \ (dR_super - dR) #Δ may be not int

		Tatompath = map(super_atomhr.atompath[i_super, j_super]) do ap
			path = TM * ap.unitpath + Δ
			if all(x -> abs(x - round(x)) ≤ 1e-3, path)
				return AtomPath(round.(Int, path), ap.hrindex)
			else
				error("cell and supercell are incompatible.")
			end
		end

		lock(lk) do
			append!(atompath[i, j], Tatompath)
		end
	end

	#Calculate some data for reindex.
	num_eachatom_orb = Vector{Int}(undef, num_atom)
	for i in 1:num_atom
		num_eachatom_orb[i] = super_atomhr.num_orb[atom_index[i][1]]
	end
	re_orb = [0; num_eachatom_orb]
	re_orb = [sum(re_orb[1:i]) for i in 1:num_atom]

	#Create cell's hr.
	hr_path = Matrix{Int}(undef, 6, 0)
	for I in CartesianIndices(atompath)
		aps = atompath[I]
		N_atompaths = map(ap -> length(ap.hrindex), aps)
		hr_path_T = Matrix{Int}(undef, 6, sum(N_atompaths))
		n = 0
		for (i, ap) in enumerate(aps)
			if N_atompaths[i] > 0
				hr_path_T[:, n+1:n+N_atompaths[i]] = [repeat([ap.unitpath; re_orb[I[1]]; re_orb[I[2]]], inner = (1, N_atompaths[i])); transpose(ap.hrindex)]
				n += N_atompaths[i]
			end
		end
		hr_path = [hr_path hr_path_T]
	end
	hr_path[4:5, :] = hr_path[4:5, :] + super_atomhr.hrorbindex[:, hr_path[6, :]]


	#orbital_index:[atom_index atom_orb_index super_atom_index super_orb_index]
	orbital_index = Matrix{Int}(undef, 0, 4)
	for i in 1:num_atom
		i_super = atom_index[i][1]
		super_orb_index = super_atomhr.orbindex[i_super]
		n = length(super_orb_index)
		orbital_index = [orbital_index; fill(i, n) collect(1:n) fill(i_super, n) super_orb_index]
	end

	return transpose(hr_path), orbital_index, atom_index
end
