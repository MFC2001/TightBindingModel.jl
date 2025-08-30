
# A symmetry (W, w) (or (S, τ)) induces a symmetry in the Brillouin zone that the
# Hamiltonian at S k is unitary equivalent to that at k, which we exploit to reduce
# computations. The relationship is
#   S = W'
#   τ = -W^-1 w
# (valid both in reduced and Cartesian coordinates). In our notation the rotation matrix
# W and translation w are such that, for each atom of type A at position a, W a + w is also
# an atom of type A.

# The full (reducible) Brillouin zone is implicitly represented by a set of (irreducible)
# kpoints (see explanation in docs). Each irreducible k-point k comes with a list of
# symmetry operations (S, τ) (containing at least the trivial operation (I, 0)), where S is
# a unitary matrix (/!\ in Cartesian but not in reduced coordinates) and τ a translation
# vector. The k-point is then used to represent implicitly the information at all the
# kpoints Sk. The relationship between the Hamiltonians is
#   H_{Sk} = U H_k U*, with
#   (Uu)(x) = u(W x + w)
# or in Fourier space
#   (Uu)(G) = e^{-i G τ} u(S^-1 G)
# In particular, we can choose the eigenvectors at Sk as u_{Sk} = U u_k

# We represent then the BZ as a set of irreducible points `kpoints`, and a set of weights
# `kweights` (summing to 1). The value of observables is given by a weighted sum over the
# irreducible kpoints, plus a symmetrization operation (which depends on the particular way
# the observable transforms under the symmetry).

# There is by decreasing cardinality
# - The group of symmetry operations of the lattice
# - The group of symmetry operations of the crystal (model.symmetries)
# - The group of symmetry operations of the crystal that preserves the BZ mesh (basis.symmetries)

# See https://juliamolsim.github.io/DFTK.jl/stable/developer/symmetries/ for details.



@doc raw"""
Return the symmetries given an atomic structure with optionally designated magnetic moments
on each of the atoms. The symmetries are determined using spglib.
"""
function symmetry_operations(lattice, atoms, positions, magnetic_moments = []; symtol = 1e-5, check_symmetry = true)
	spin_polarization = determine_spin_polarization(magnetic_moments)
	dimension = count(!iszero, eachcol(lattice))
	if isempty(atoms) || dimension != 3
		# spglib doesn't support these cases, so we default to no symmetries
		return [one(SymOp)]
	end

	if spin_polarization == :full
		@warn("Symmetry detection not yet supported in full spin polarization. " *
			  "Returning no symmetries")
		return [one(SymOp)]
	end

	atom_groups = [findall(Ref(pot) .== atoms) for pot in Set(atoms)]
	cell = spglib_cell(lattice, atom_groups, positions, magnetic_moments)
	(Ws, ws) = try
		if spin_polarization == :none
			Spglib.get_symmetry(cell, symtol)
		elseif spin_polarization == :collinear
			Spglib.get_symmetry_with_collinear_spin(cell, symtol)
		end
	catch e
		if e isa Spglib.SpglibError
			msg = ("spglib failed to get the symmetries. Check your lattice, use a " *
				   "uniform BZ mesh or disable symmetries. Spglib reported : " * e.msg)
			throw(Spglib.SpglibError(msg))
		else
			rethrow()
		end
	end

	symmetries = [SymOp(W, w) for (W, w) in zip(Ws, ws)]
	if check_symmetry
		_check_symmetries(symmetries, lattice, atom_groups, positions; symtol)
	end
	return symmetries
end

"""
Return the Symmetry operations given a `hall_number`.

This function allows to directly access to the space group operations in the
`spglib` database. To specify the space group type with a specific choice,
`hall_number` is used.

The definition of `hall_number` is found at
[Space group type](https://spglib.readthedocs.io/en/latest/dataset.html#dataset-spg-get-dataset-spacegroup-type).
"""
function symmetry_operations(hall_number::Integer)
	(Ws, ws) = Spglib.get_symmetry_from_database(hall_number)
	return [SymOp(W, w) for (W, w) in zip(Ws, ws)]
end
