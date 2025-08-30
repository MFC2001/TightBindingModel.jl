export TwoLayerPOSCAR2HR
function TwoLayerPOSCAR2HR(
	unithr::HR,
	unitorbital::ORBITAL,
	unitposcar::POSCAR,
	poscar::POSCAR,
	zinterface::Real,
	interlayerU::Function;
	zunit = "Cartesian",
	finduc::AbstractString = "auto",
	uc_atom_orb::Real = 0.5,
	interlayerUeps::Real = 1e-6,
	maxdist::Real = 0,
)

	layerposcar = unstackPOSCAR(poscar, zinterface; unit = zunit)

	nlayer = 2
	layerhr = Vector{HR}(undef, nlayer)
	layerorbital = Vector{ORBITAL}(undef, nlayer)
	for i in eachindex(layerposcar)
		(layerhr[i], layerorbital[i]) = POSCAR2HR(unithr, unitorbital, unitposcar, layerposcar[i]; finduc, outhr = "value", outorb = "value", uc_atom_orb)
	end


	reindexHR!(layerhr[2], layerorbital[1].num)


	#make layerorbital's index for building interlayer HR.
	layerorbitalindex = Vector{Vector{Int}}(undef, nlayer)
	layerorbitalindex[1] = (1:layerorbital[1].num)
	layerorbitalindex[2] = (1:layerorbital[2].num) .+ layerorbital[1].num

	interlayerhr = ORBITAL2HR(poscar, layerorbital[1], layerorbitalindex[1], layerorbital[2], layerorbitalindex[2], interlayerU; heps = interlayerUeps, maxdist)

	hr = integrate(layerhr..., interlayerhr)
	orbital = integrate(layerorbital...)

	return hr, orbital
end
