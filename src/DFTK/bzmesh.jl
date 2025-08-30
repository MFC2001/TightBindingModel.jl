
# import Brillouin.KPaths: KPathInterpolant




@doc raw"""
Build a [`MonkhorstPack`](@ref) grid to ensure kpoints are at most this `spacing`
apart (in inverse Bohrs). A reasonable spacing is `0.13` inverse Bohrs
(around ``2π * 0.04 \, \text{Å}^{-1}``). The `kshift` keyword argument allows
to specify an explicit shift for all ``k``-points.
"""
# function kgrid_from_maximal_spacing(system::AbstractSystem, spacing; kshift=[0, 0, 0])
#     pseudopotentials = fill(nothing, length(system))
#     kgrid_from_maximal_spacing(parse_system(system, pseudopotentials).lattice, spacing; kshift)
# end
function kgrid_from_maximal_spacing(lattice::AbstractMatrix, spacing; kshift = [0, 0, 0])
	lattice       = austrip.(lattice)
	spacing       = austrip(spacing)
	recip_lattice = compute_recip_lattice(lattice)
	@assert spacing > 0
	isinf(spacing) && return MonkhorstPack([1, 1, 1], kshift)

	kgrid = [max(1, ceil(Int, norm(vec) / spacing)) for vec in eachcol(recip_lattice)]
	MonkhorstPack(kgrid, kshift)
end
# function kgrid_from_maximal_spacing(model::Model, spacing; kwargs...)
#     kgrid_from_maximal_spacing(model.lattice, spacing; kwargs...)
# end

@doc raw"""
Selects a [`MonkhorstPack`](@ref) grid size which ensures that at least a
`n_kpoints` total number of ``k``-points are used. The distribution of
``k``-points amongst coordinate directions is as uniformly as possible, trying to
achieve an identical minimal spacing in all directions.
"""
# function kgrid_from_minimal_n_kpoints(system::AbstractSystem, n_kpoints::Integer; kshift=[0, 0, 0])
#     pseudopotentials = fill(nothing, length(system))
#     kgrid_from_minimal_n_kpoints(parse_system(system, pseudopotentials).lattice, n_kpoints; kshift)
# end
function kgrid_from_minimal_n_kpoints(lattice, n_kpoints::Integer; kshift = [0, 0, 0])
	lattice = austrip.(lattice)
	n_dim   = count(!iszero, eachcol(lattice))
	@assert n_kpoints > 0
	n_kpoints == 1 && return MonkhorstPack([1, 1, 1], kshift)
	n_dim == 1 && return MonkhorstPack([n_kpoints, 1, 1], kshift)

	# Compute truncated reciprocal lattice
	recip_lattice_nD = 2π * inv(lattice[1:n_dim, 1:n_dim]')
	n_kpt_per_dim = n_kpoints^(1 / n_dim)

	# Start from a cubic lattice. If it is one, we are done. Otherwise the resulting
	# spacings in each dimension bracket the ideal k-point spacing.
	spacing_per_dim = [norm(vec) / n_kpt_per_dim for vec in eachcol(recip_lattice_nD)]
	min_spacing, max_spacing = extrema(spacing_per_dim)
	if min_spacing ≈ max_spacing
		return kgrid_from_maximal_spacing(lattice, min_spacing)
	else
		number_of_kpoints(spacing) = prod(vec -> norm(vec) / spacing, eachcol(recip_lattice_nD))
		@assert number_of_kpoints(min_spacing) + 0.05 ≥ n_kpoints
		@assert number_of_kpoints(max_spacing) - 0.05 ≤ n_kpoints

		# TODO This is not optimal and sometimes finds too large grids
		spacing = Roots.find_zero(sp -> number_of_kpoints(sp) - n_kpoints,
			(min_spacing, max_spacing), Roots.Bisection(),
			xatol = 1e-4, atol = 0, rtol = 0)

		# Sanity check: Sometimes root finding is just across the edge towards
		# a larger number of k-points than needed. This attempts a slightly larger spacing.
		kgrid_larger = kgrid_from_maximal_spacing(lattice, spacing + 1e-4)
		if length(kgrid_larger) ≥ n_kpoints
			return kgrid_larger
		else
			return kgrid_from_maximal_spacing(lattice, spacing)
		end
	end
end
# function kgrid_from_minimal_n_kpoints(model::Model, n_kpoints::Integer; kwargs...)
#     kgrid_from_minimal_n_kpoints(model.lattice, n_kpoints; kwargs...)
# end


