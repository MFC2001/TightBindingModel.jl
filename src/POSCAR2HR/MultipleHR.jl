
function ExpandHR(
	unithr::AbstractVector{<:HR},
	unitorbital::ORBITAL,
	unitcell::Cell,
	cell::Cell; #larger than unitcell.
	finduc::AbstractString = "auto",
	supercell_path::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
	ucazimuth::Tuple{<:Real, <:Real} = (0, 0),
	outhr::AbstractString = "value",
	outorb::AbstractString = "value",
	outucorientation::AbstractString = "N",
)

	if outucorientation[1] ∈ ['N', 'n']
		(hr_index, orbital_index) = ExpandHR(unithr[1], unitorbital, unitcell, cell; finduc, supercell_path, ucazimuth, outhr = "index", outorb = "index", outucorientation)
        
        for i in eachindex(unithr)
            
        end
	end



	ExpandHR(unithr, unitorbital, unitcell, cell; finduc, supercell_path, ucazimuth, outhr = "index", outorb = "index", outucorientation)


end
