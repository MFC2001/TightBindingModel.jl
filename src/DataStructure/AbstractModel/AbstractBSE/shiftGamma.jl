function _BSE_shiftΓ(period)
	if count(period) == 0
		q = ReducedCoordinates(0, 0, 0)
	else
		I = findfirst(period)
		q = [0.0, 0.0, 0.0]
		q[I] = 1e-8
		q = ReducedCoordinates(q)
	end
	return q
end
