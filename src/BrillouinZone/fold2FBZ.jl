export fold2FBZ

function fold2FBZ(Rlattice::ReciprocalLattice, kpoint::AbstractVector{<:AbstractVector})

	equalΓ = [
		[-1, -1, -1],
		[-1, -1, 0],
		[-1, -1, 1],
		[-1, 0, -1],
		[-1, 0, 0],
		[-1, 0, 1],
		[-1, 1, -1],
		[-1, 1, 0],
		[-1, 1, 1],
		[0, -1, -1],
		[0, -1, 0],
		[0, -1, 1],
		[0, 0, -1],
		[0, 0, 1],
		[0, 1, -1],
		[0, 1, 0],
		[0, 1, 1],
		[1, -1, -1],
		[1, -1, 0],
		[1, -1, 1],
		[1, 0, -1],
		[1, 0, 0],
		[1, 0, 1],
		[1, 1, -1],
		[1, 1, 0],
		[1, 1, 1],
	]


	foldedkpoint = deepcopy(kpoint)
	for kindex in eachindex(foldedkpoint)

		kpoint = foldedkpoint[kindex]

		if _is_approx_integer(kpoint; atol = 1e-4)
			foldedkpoint[kindex] = kpoint - round.(kpoint)
			continue
		end

		kpoint_coord = Rlattice * kpoint
	
		notinFBZ = true
		while notinFBZ
			Δ2Γ = norm(kpoint_coord)
			notinFBZ = false
			for eΓ in equalΓ
				T = kpoint_coord - Rlattice * eΓ
				if norm(T) < Δ2Γ
					kpoint_coord = T
					kpoint = kpoint - eΓ
					notinFBZ = true
					break
				end
			end
		end

		foldedkpoint[kindex] = kpoint
	end

	return foldedkpoint
end

function fold2FBZ(Rlattice::ReciprocalLattice, kpoint::AbstractVector{<:Real})

	if _is_approx_integer(kpoint; atol = 1e-4)
		return kpoint - round.(kpoint)
	end

	equalΓ = [
		[-1, -1, -1],
		[-1, -1, 0],
		[-1, -1, 1],
		[-1, 0, -1],
		[-1, 0, 0],
		[-1, 0, 1],
		[-1, 1, -1],
		[-1, 1, 0],
		[-1, 1, 1],
		[0, -1, -1],
		[0, -1, 0],
		[0, -1, 1],
		[0, 0, -1],
		[0, 0, 1],
		[0, 1, -1],
		[0, 1, 0],
		[0, 1, 1],
		[1, -1, -1],
		[1, -1, 0],
		[1, -1, 1],
		[1, 0, -1],
		[1, 0, 0],
		[1, 0, 1],
		[1, 1, -1],
		[1, 1, 0],
		[1, 1, 1],
	]

	kpoint_coord = Rlattice * kpoint

	notinFBZ = true
	while notinFBZ
		Δ2Γ = norm(kpoint_coord)
		notinFBZ = false
		for eΓ in equalΓ
			T = kpoint_coord - Rlattice * eΓ
			if norm(T) < Δ2Γ
				kpoint_coord = T
				kpoint = kpoint - eΓ
				notinFBZ = true
				break
			end
		end
	end

	return kpoint
end
