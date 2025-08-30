@doc raw"""
Apply various standardisations to a lattice and a list of atoms. It uses spglib to detect
symmetries (within `tol_symmetry`), then cleans up the lattice according to the symmetries
(unless `correct_symmetry` is `false`) and returns the resulting standard lattice
and atoms. If `primitive` is `true` (default) the primitive unit cell is returned, else
the conventional unit cell is returned.
"""
function standardize_atoms(lattice, atoms, positions, magnetic_moments = []; kwargs...)
	@assert length(atoms) == length(positions)
	@assert isempty(magnetic_moments) || (length(atoms) == length(magnetic_moments))
	atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]
	ret = spglib_standardize_cell(lattice, atom_groups, positions, magnetic_moments; kwargs...)
	(; ret.lattice, atoms, ret.positions, ret.magnetic_moments)
end

"""
Returns crystallographic conventional cell according to the International Table of
Crystallography Vol A (ITA) in case `primitive=false`. If `primitive=true`
the primitive lattice is returned in the convention of the reference work of
Cracknell, Davies, Miller, and Love (CDML). Of note this has minor differences to
the primitive setting choice made in the ITA.
"""
function spglib_standardize_cell(lattice::AbstractArray{T}, atom_groups, positions, magnetic_moments = [];
	correct_symmetry = true, primitive = false, symtol = 1e-5) where {T}
	# TODO For time-reversal symmetry see the discussion in PR 496.
	#      https://github.com/JuliaMolSim/DFTK.jl/pull/496/files#r725203554
	#      Essentially this does not influence the standardisation,
	#      but it only influences the kpath.
	cell = spglib_cell(lattice, atom_groups, positions, magnetic_moments)
	std_cell = Spglib.standardize_cell(cell, tol_symmetry = symtol; to_primitive = primitive,
		no_idealize = !correct_symmetry)

	lattice = Matrix{T}(std_cell.lattice)
	positions = Vec3{T}.(std_cell.positions)
	magnetic_moments = normalize_magnetic_moment.(std_cell.magmoms)
	(; lattice, atom_groups, positions, magnetic_moments)
end
