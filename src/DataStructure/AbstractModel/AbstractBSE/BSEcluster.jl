struct BSEcluster <: AbstractBSE

end




function BSE(unitTB::AbstractTightBindModel{HT}, supercell::Cell,
	kgrid::MonkhorstPack, vindex::Vector{<:Integer}, cindex::Vector{<:Integer};
	finduc::AbstractString = "auto",
	scissor::Real = 0,
	αrcut::Real = 4.5,
	δ::Real = 1e-6,
	U_screened::Union{HR, Nothing} = nothing,
	J_screened::Union{HR, Nothing} = nothing,
	A_screened::Union{HR, Nothing} = nothing,
	U_unscreened::Union{HR, Nothing} = nothing,
	J_unscreened::Union{HR, Nothing} = nothing,
	A_unscreened::Union{HR, Nothing} = nothing,
) where (HT <: Number)

	unitcell = Cell(unitTB)
	unithr = HR(unitTB)
	unitorbital = ORBITAL(TB)
	(unitcell, supercell) = POSCAR2HRs.checkCell(unitcell, supercell)

	if finduc ∉ ["auto", "translate", "custom"]
		error("Wrong keyword finduc.")
	end

	EpsPara = Dict(
		#Used by classifyorbital, permissible error between orbital center location and atom location.
		# "uc atom_orb" => 0.5,
		#Used by findunitcell!, permissible error between atom locations in operated unitcell and cell's supercell.
		"atomeps" => 0.1,
		#Used by findunitcell!, sum square of permissible error between atom locations in operated unitcell and cell's supercell.
		"sumeps" => 0.2,
		#Used by classifyPOSCAR, permissible error between atom locations in unitcell's supercell and cell.
		"sc_uc atom" => 0.2,
		#Used by atompath_sum in poscar2hr, permissible error when search cell's atom transition path in unithr.path.
		"sc_uc path" => 0.2,
		#Used by maxhrpath in hrsplit, the added value of the farthest transition lattice distance in unithr.path.
		"unithr path" => 0.1,
	)

	#This function will use name to judge whether atom is the same, so need cell.name is corresponding to unitcell.name, recommend using atom name.
	if finduc == "custom"
	else
		(unitcell, ucazimuth) = POSCAR2HRs.findunitcell(cell, unitcell; aimazimuth = ucazimuth, mode = finduc, atomeps = EpsPara["atomeps"], sumeps = EpsPara["atomeps"])
	end
	#Preprocessing unitcell HR.
	POSCAR2HRs.namelist!(unitcell)

	unit_atomhr = POSCAR2HRs.hrsplit(unithr, unitcell, unitorbital.belonging, EpsPara["unithr path"])
	(hr_path, orbital_index, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atomhr, EpsPara)

	superhr = HR(hr_path[:, 1:5], unithr.value[hr_path[:, 6]]; hrsort = "Y")

	orbital_dlocation_frac = POSCAR2HRs.ORBITAL_fracdlocation(unitorbital, unitcell.lattice)
	orbital_name = cell.name[orbital_index[:, 1]]
	orbital_location = map((ai, uoi) -> cell.location[ai] + unitcell.lattice * orbital_dlocation_frac[uoi], orbital_index[:, 1], orbital_index[:, 4])
	superorbital = ORBITAL(orbital_location; name = orbital_name, index = orbital_index[:, 4],
		atom_location = cell.location, atom_name = cell.name, belonging = orbital_index[:, 1])

	superTB = TightBindModel(superhr, supercell, superorbital)

	norb = length(orbital_location)
	ZERO = zeros(Int, norb, norb)

	if isnothing(U_screened)
		println("Direct part of W in K^d is zero!")
		U_screened = (k) -> ZERO
	else
		unit_atom_U_screened = POSCAR2HRs.hrsplit(U_screened, unitcell, unitorbital.belonging, EpsPara["unithr path"])
		(hr_path, _, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atom_U_screened, EpsPara)
		U_screened = HR(hr_path[:, 1:5], U_screened.value[hr_path[:, 6]]; hrsort = "Y")
		U_screened = VR2VK(U_screened)
	end
	if isnothing(J_screened)
		J_screened = (k) -> ZERO
	elseif J_screened isa HR
		unit_atom_J_screened = POSCAR2HRs.hrsplit(J_screened, unitcell, unitorbital.belonging, EpsPara["unithr path"])
		(hr_path, _, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atom_J_screened, EpsPara)
		J_screened = HR(hr_path[:, 1:5], J_screened.value[hr_path[:, 6]]; hrsort = "Y")
		J_screened = VR2VK(J_screened)
	end
	if isnothing(A_screened)
		A_screened = (k) -> ZERO
	elseif A_screened isa HR
		unit_atom_A_screened = POSCAR2HRs.hrsplit(A_screened, unitcell, unitorbital.belonging, EpsPara["unithr path"])
		(hr_path, _, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atom_A_screened, EpsPara)
		A_screened = HR(hr_path[:, 1:5], A_screened.value[hr_path[:, 6]]; hrsort = "Y")
		A_screened = VR2VK(A_screened)
	end
	if isnothing(J_unscreened)
		J_unscreened = (k) -> ZERO
	elseif J_unscreened isa HR
		unit_atom_J_unscreened = POSCAR2HRs.hrsplit(J_unscreened, unitcell, unitorbital.belonging, EpsPara["unithr path"])
		(hr_path, _, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atom_J_unscreened, EpsPara)
		J_unscreened = HR(hr_path[:, 1:5], J_unscreened.value[hr_path[:, 6]]; hrsort = "Y")
		J_unscreened = VR2VK(J_unscreened)
	end
	if isnothing(A_unscreened)
		A_unscreened = (k) -> ZERO
	elseif A_unscreened isa HR
		unit_atom_A_unscreened = POSCAR2HRs.hrsplit(A_unscreened, unitcell, unitorbital.belonging, EpsPara["unithr path"])
		(hr_path, _, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atom_A_unscreened, EpsPara)
		A_unscreened = HR(hr_path[:, 1:5], A_unscreened.value[hr_path[:, 6]]; hrsort = "Y")
		A_unscreened = VR2VK(A_unscreened)
	end
	if isnothing(U_unscreened)
		U_unscreened = (k) -> ZERO
		println("Direct part of V in K^x is zero!")
	elseif U_unscreened isa HR
		unit_atom_U_unscreened = POSCAR2HRs.hrsplit(U_unscreened, unitcell, unitorbital.belonging, EpsPara["unithr path"])
		(hr_path, _, _) = POSCAR2HRs.poscar2hr(cell, unitcell, unit_atom_U_unscreened, EpsPara)
		U_unscreened = HR(hr_path[:, 1:5], U_unscreened.value[hr_path[:, 6]]; hrsort = "Y")
		U_unscreened = _BSE_exchange(superTB, U_unscreened, αrcut, δ)
	end


	nv = length(vindex)
	nc = length(cindex)
	Nk = length(kgrid)
	N = nv * nc * Nk
	vckindex = Vector{Tuple{Int, Int, Int}}(undef, N)
	n = 0
	for k in 1:Nk, c in 1:nc, v in 1:nv
		n += 1
		vckindex[n] = (v, c, k)
	end


	return BSE{HT, typeof(superTB)}(superTB, kgrid, vindex, cindex, vckindex, scissor,
		U_screened, J_screened, A_screened, U_unscreened, J_unscreened, A_unscreened)
end


