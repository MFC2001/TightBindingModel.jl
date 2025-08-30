"""
	ExpandHR(uchrfile::String, ucorbitalfile, ucposcarfile, scposcarfile;
		heps::Real = 0,	fermienergy::Real = 0, ucperiodicity = Bool[1, 1, 1],	scperiodicity = Bool[1, 1, 1], 
		finduc = "auto", supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing)
"""
function ExpandHR(
	uchrfile::AbstractString,
	ucorbitalfile::AbstractString,
	ucposcarfile::AbstractString,
	scposcarfile::AbstractString;
	heps::Real = 0,
	fermienergy::Real = 0,
	ucperiodicity = [1, 1, 1],
	scperiodicity = [1, 1, 1],
	finduc::AbstractString = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
)

	uchr = ReadHR(uchrfile, heps; readimag = 'Y', μ = fermienergy)
	ucorbital = ReadORBITAL(ucorbitalfile)
	uccell = ReadPOSCAR(ucposcarfile; period = ucperiodicity)

	sccell = ReadPOSCAR(scposcarfile; period = scperiodicity)

	(schr, scorbital) = ExpandHR(uchr, ucorbital, uccell, sccell; finduc, supercell_path)

	scfolder = dirname(scposcarfile)
	WriteHR(schr, joinpath(scfolder, "hr.dat"))
	WriteORBITAL(scorbital, joinpath(scfolder, "orbital.xyz"))

	return nothing
end
"""
	ExpandHR(unitTB::AbstractTightBindModel, supercell::Cell;
		finduc::AbstractString = "auto", supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing, ucazimuth::Tuple{<:Real, <:Real} = (0, 0), 
		outucorientation::AbstractString = "N")
"""
function ExpandHR(
	unitTB::AbstractTightBindModel,
	supercell::Cell;
	finduc::AbstractString = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
	ucazimuth::Tuple{<:Real, <:Real} = (0, 0),
	outucorientation::AbstractString = "N",
)

	unitcell = Cell(unitTB)
	unithr = HR(unitTB)
	unitorbital = ORBITAL(unitTB)

	return ExpandHR(unithr, unitorbital, unitcell, supercell;
		finduc, supercell_path, ucazimuth, outhr = "value", outorb = "value", outucorientation)
end

"""
	ContractHR(schrfile::String, scorbitalfile, scposcarfile, ucposcarfile;
		heps::Real = 0,	fermienergy::Real = 0, ucperiodicity = Bool[1, 1, 1],	scperiodicity = Bool[1, 1, 1], 
		finduc = "auto", supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing)
"""
function ContractHR(
	schrfile::AbstractString,
	scorbitalfile::AbstractString,
	scposcarfile::AbstractString,
	ucposcarfile::AbstractString;
	heps::Real = 0,
	fermienergy::Real = 0,
	ucperiodicity = [1, 1, 1],
	scperiodicity = [1, 1, 1],
	finduc::AbstractString = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
)

	schr = ReadHR(schrfile, heps; readimag = 'Y', μ = fermienergy)
	scorbital = ReadORBITAL(scorbitalfile)
	sccell = ReadPOSCAR(scposcarfile; period = scperiodicity)

	uccell = ReadPOSCAR(ucposcarfile; period = ucperiodicity)

	(uchr, ucorbital) = ContractHR(schr, scorbital, sccell, uccell; finduc, supercell_path)

	ucfolder = dirname(ucposcarfile)
	WriteHR(uchr, joinpath(ucfolder, "hr.dat"))
	WriteORBITAL(ucorbital, joinpath(ucfolder, "orbital.xyz"))

	return nothing
end
"""
	ContractHR(superTB::AbstractTightBindModel, unitcell::Cell;
		finduc::AbstractString = "auto", supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing, ucazimuth::Tuple{<:Real, <:Real} = (0, 0), 
		outucorientation::AbstractString = "N")
"""
function ContractHR(
	superTB::AbstractTightBindModel,
	unitcell::Cell;
	finduc::AbstractString = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
	ucazimuth::Tuple{<:Real, <:Real} = (0, 0),
	outucorientation::AbstractString = "N",
)

	supercell = Cell(superTB)
	superhr = HR(superTB)
	superorbital = ORBITAL(superTB)

	return ContractHR(superhr, superorbital, supercell, unitcell;
		finduc, supercell_path, ucazimuth, outhr = "value", outorb = "value", outucorientation)
end
