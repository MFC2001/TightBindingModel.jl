
normalize_magnetic_moment(::Nothing)::Vec3{Float64}          = (0, 0, 0)
normalize_magnetic_moment(mm::Number)::Vec3{Float64}         = (0, 0, mm)
normalize_magnetic_moment(mm::AbstractVector)::Vec3{Float64} = mm

"""
`:none` if no element has a magnetic moment, else `:collinear` or `:full`.
"""
function determine_spin_polarization(magnetic_moments)
	isempty(magnetic_moments) && return :none
	all_magmoms = normalize_magnetic_moment.(magnetic_moments)
	all(iszero, all_magmoms) && return :none
	all(iszero(magmom[1:2]) for magmom in all_magmoms) && return :collinear

	return :full
end