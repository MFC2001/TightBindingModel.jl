# Routines for interaction with spglib
#
# Note that these routines should be generalised to be compatible with AtomsBase structures
# and not live in DFTK, but we keep them for now.

function spglib_cell(lattice, atom_groups, positions, magnetic_moments)
	@assert !isempty(atom_groups)  # Otherwise spglib cannot work properly
	magnetic_moments = normalize_magnetic_moment.(magnetic_moments)

	spg_atoms     = Int[]
	spg_magmoms   = Float64[]
	spg_positions = Vector{Float64}[]

	arbitrary_spin = false
	for (igroup, indices) in enumerate(atom_groups), iatom in indices
		push!(spg_atoms, igroup)
		push!(spg_positions, positions[iatom])

		if isempty(magnetic_moments)
			magmom = zeros(3)
		else
			magmom = magnetic_moments[iatom]
			!iszero(magmom[1:2]) && (arbitrary_spin = true)
		end
		push!(spg_magmoms, magmom[3])
	end
	@assert !arbitrary_spin
	Spglib.SpglibCell(lattice, spg_positions, spg_atoms, spg_magmoms)
end
