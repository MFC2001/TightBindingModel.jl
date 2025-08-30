function _check_symmetries(symmetries::AbstractVector{<:SymOp}, lattice, atom_groups, positions; symtol = 1e-5)
	# Check (W, w) maps atoms to equivalent atoms in the lattice
	for symop in symmetries
		W, w = symop.W, symop.w

		# Check (A W A^{-1}) is orthogonal
		Wcart = lattice * W / lattice
		if maximum(abs, Wcart'Wcart - I) > symtol
			error("Issue in symmetry determination: Non-orthogonal rotation matrix.")
		end

		for group in atom_groups
			group_positions = positions[group]
			for coord in group_positions
				# If all elements of a difference in diffs is integer, then
				# W * coord + w and pos are equivalent lattice positions
				if !any(c -> _is_approx_integer(W * coord + w - c; atol = symtol), group_positions)
					error("Issue in symmetry determination: Cannot map the atom at position " *
						  "$coord to another atom of the same element under the symmetry " *
						  "operation (W, w):\n($W, $w)")
				end
			end
		end
	end
end

function _check_kpoint_reduction(symmetries::AbstractVector{<:SymOp}, kgrid::MonkhorstPack, k_all_reducible, kirreds, grid_address)
	for (iks_reducible, k) in zip(k_all_reducible, kirreds), ikred in iks_reducible
		kred = (kgrid.kshift .+ grid_address[ikred]) .// kgrid.kgrid_size
		found_mapping = any(symmetries) do symop
			# If the difference between kred and W' * k == W^{-1} * k
			# is only integer in fractional reciprocal-space coordinates, then
			# kred and S' * k are equivalent k-points
			all(isinteger, kred - (symop.S * k))
		end
		if !found_mapping
			error("The reducible k-point $kred could not be generated from " *
				  "the irreducible kpoints. This points to a bug in spglib.")
		end
	end
end
