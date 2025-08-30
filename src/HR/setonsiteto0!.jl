function setonsiteto0!(Uhr::HR{T}) where {T}

	I = findall(x -> iszero(x[1:3]) && isequal(x[4], x[5]), eachrow(Uhr.path))
	Uhr.value[I] .= 0

	for i in eachindex(Uhr.orbindex)
		for (ii, hop) in enumerate(Uhr.hop[i, i])
			if iszero(hop.R)
				Uhr.hop[i, i][ii] = Hopping{T}(hop.i, hop.j, hop.R, 0)
			end
		end
	end

	return Uhr
end
