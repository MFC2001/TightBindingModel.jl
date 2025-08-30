
"""
Filter out the symmetry operations that don't respect the symmetries of the discrete BZ grid.
"""
function symmetries_preserving_kgrid(symmetries, kcoords)
	kcoords_normalized = normalize_kpoint_coordinate.(kcoords)
	function preserves_grid(symop)
		all(_is_approx_in(normalize_kpoint_coordinate(symop.S * k), kcoords_normalized) for k in kcoords_normalized)
	end
	return filter(preserves_grid, symmetries)
end
function symmetries_preserving_kgrid(symmetries, kgrid::ExplicitKpoints; symtol = 1e-5)
	# First apply symmetries as the provided k-points can be arbitrary
	# (e.g. only along a line or similar)
	all_kcoords = unfold_kcoords(kgrid.kcoords, symmetries; symtol)
	return symmetries_preserving_kgrid(symmetries, all_kcoords)
end
function symmetries_preserving_kgrid(symmetries, kgrid::MonkhorstPack)
	if all(isone, kgrid.kgrid_size)
		# TODO Keeping this special casing from version of the code before refactor
		return [one(SymOp)]
	else
		# if k' = Rk
		# then
		#    R' = diag(kgrid) R diag(kgrid)^-1
		# should be integer where

		# TODO This can certainly be improved by knowing this is an MP grid,
		#      see symmetries_preserving_rgrid below for ideas
		return symmetries_preserving_kgrid(symmetries, reducible_kcoords(kgrid).kcoords)
	end
end

# Approximate in; can be performance-critical, so we optimize in case of rationals
_is_approx_in(x::AbstractArray{<:Rational}, X)  = any(isequal(x), X)
_is_approx_in(x::AbstractArray{T}, X) where {T} = any(y -> isapprox(x, y; atol = sqrt(eps(T))), X)




"""
Filter out the symmetry operations that don't respect the symmetries of the discrete real-space grid
"""
function symmetries_preserving_rgrid(symmetries, fft_size; symtol = 1e-5)
	is_in_grid(r) =
		all(zip(r, fft_size)) do (ri, size)
			abs(ri * size - round(ri * size)) / size ≤ symtol
		end

	onehot3(i) = (x = zeros(Bool, 3); x[i] = true; Vec3(x))
	function preserves_grid(symop)
		all(is_in_grid(symop.W * onehot3(i) .// fft_size[i] + symop.w) for i ∈ 1:3)
	end

	filter(preserves_grid, symmetries)
end
