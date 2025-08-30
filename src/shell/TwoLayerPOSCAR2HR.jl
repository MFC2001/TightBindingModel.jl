export TwoLayerPOSCAR2HR
function TwoLayerPOSCAR2HR(
	uchrfile::AbstractString,
	ucorbitalfile::AbstractString,
	ucposcarfile::AbstractString,
	scposcarfile::AbstractString,
	interlayerU::Function;
	heps::Real = 1e-6,
	fermienergy::Real = 0,
	ucperiodicity = ["p", "p", "p"],
	scperiodicity = ["p", "p", "p"],
	finduc::AbstractString = "auto",
	zinterface::Real = 0,
	maxdist::Real = 0,
)

	uchr = ReadHR(uchrfile, heps; fermienergy)
	ucorbital = ReadORBITAL(ucorbitalfile)
	ucposcar = ReadPOSCAR(ucposcarfile; periodicity = ucperiodicity)


	scposcar = ReadPOSCAR(scposcarfile; periodicity = scperiodicity)

	(schr, scorbital) = TwoLayerPOSCAR2HR(uchr, ucorbital, ucposcar, scposcar, zinterface, interlayerU; zunit = "C", finduc, uc_atom_orb = 0.5, interlayerUeps = heps, maxdist)


	scfolder = dirname(scposcarfile)
	WriteHR(schr, joinpath(scfolder, "hr.dat"))
	WriteORBITAL(scorbital, joinpath(scfolder, "orbital.xyz"))

	return nothing
end
