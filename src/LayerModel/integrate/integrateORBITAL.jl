export integrate
"""
	integrateORBITAL(orbitals::ORBITAL...)::ORBITAL

P.S. This function will try to create a new orbital with the same type of orbitals[1], so make sure that all hrs satisfy the condition.
"""
function integrate(orbitals::ORBITAL...)::ORBITAL

	orblocation = Vector{eltype(orbitals[1].location)}(undef, 0)
	atomlocation = Vector{eltype(orbitals[1].atomlocation)}(undef, 0)
	atomname = Vector{String}(undef, 0)

	for orbital in orbitals
		append!(orblocation, orbital.location)
		append!(atomlocation, orbital.atomlocation)
		append!(atomname, orbital.atomname)
	end

	return ORBITAL(orblocation; atomlocation, atomname)
end
