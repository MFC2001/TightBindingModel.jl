export WilsonLoop
function WilsonLoop(TB::AbstractTightBindModel, bandindex::Integer = 1; Nkx = 51, Nky = 51) where {T}

	start = -floor.(Int, ([Nkx, Nky] .- 1) .// 2)
	stop = ceil.(Int, ([Nkx, Nky] .- 1) .// 2)

	#W(ky)
	ky = start[2]-1:stop[2]
	θ_ky = Vector{Float64}(undef, length(ky))
	Threads.@threads for kyi in eachindex(ky)
		aimkdirects = [Vec3([kx // Nkx, ky[kyi] // Nky, 0]) for kx in start[1]:stop[1]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nkx-1
			uu = uband[i+1].vectors[:, bandindex] ⋅ uband[i].vectors[:, bandindex]
			W *= uu
		end
		uu = sum(i -> cis(2π * ([1, 0, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, bandindex]) * uband[end].vectors[i, bandindex], eachindex(TB.orb_location))
		W = uu * W

		θ_ky[kyi] = angle(W)
	end
	θ_ky = smooth_AngleSingular(θ_ky, π)
	wcc_ky = θ_ky ./ 2π

	#W(kx)
	kx = start[1]-1:stop[1]
	θ_kx = Vector{Float64}(undef, length(kx))
	Threads.@threads for kxi in eachindex(kx)
		aimkdirects = [Vec3([kx[kxi] // Nkx, ky // Nky, 0]) for ky in start[2]:stop[2]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nky-1
			uu = uband[i+1].vectors[:, bandindex] ⋅ uband[i].vectors[:, bandindex]
			W *= uu
		end
		uu = sum(i -> cis(2π * ([0, 1, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, bandindex]) * uband[end].vectors[i, bandindex], eachindex(TB.orb_location))
		W = uu * W

		θ_kx[kxi] = angle(W)
	end
	θ_kx = smooth_AngleSingular(θ_kx, π)
	wcc_kx = θ_kx ./ 2π

	return wcc_kx, wcc_ky
end

function WilsonLoop(TB::AbstractTightBindModel, bandindex::AbstractVector{<:Integer} = [1, 2]; Nkx = 51, Nky = 51) where {T}

	start = -floor.(Int, ([Nkx, Nky] .- 1) .// 2)
	stop = ceil.(Int, ([Nkx, Nky] .- 1) .// 2)

	nband = length(bandindex)

	#W(ky)
	ky = start[2]-1:stop[2]
	θ_ky = Matrix{Float64}(undef, nband, length(ky))
	Threads.@threads for kyi in eachindex(ky)
		aimkdirects = [Vec3([kx // Nkx, ky[kyi] // Nky, 0]) for kx in start[1]:stop[1]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nkx-1
			uu = [uband[i+1].vectors[:, m] ⋅ uband[i].vectors[:, n] for m in bandindex, n in bandindex]
			F = svd(uu)
			uu = F.U * F.Vt
			W = uu * W
		end
		uu = [sum(i -> cis(2π * ([1, 0, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, m]) * uband[end].vectors[i, n], eachindex(TB.orb_location)) for m in bandindex, n in bandindex]
		W = uu * W

		θ_ky[:, kyi] = eigvals!(Hermitian(W))
	end
	wcc_ky = θ_ky ./ 2π


	#W(kx)
	kx = start[1]-1:stop[1]
	θ_kx = Matrix{Float64}(undef, nband, length(kx))
	Threads.@threads for kxi in eachindex(kx)
		aimkdirects = [Vec3([kx[kxi] // Nkx, ky // Nky, 0]) for ky in start[2]:stop[2]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nky-1
			uu = [uband[i+1].vectors[:, m] ⋅ uband[i].vectors[:, n] for m in bandindex, n in bandindex]
			F = svd(uu)
			uu = F.U * F.Vt
			W = uu * W
		end
		uu = [sum(i -> cis(2π * ([0, 1, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, m]) * uband[end].vectors[i, n], eachindex(TB.orb_location)) for m in bandindex, n in bandindex]
		W = uu * W

		θ_kx[kxi] = eigvals!(Hermitian(W))
	end
	wcc_kx = θ_kx ./ 2π

	return wcc_kx, wcc_ky
end

# function _wilsonloop_uband_to_W()
# 	uband = BAND(kdirects, TB, TB.orb_location; vector = true)

# 	Nk = length(kdirects)

# 	if bandindex isa Integer
# 		W = 1
# 		for i in 1:Nk-1
# 			uu = uband[i+1].vectors[:, bandindex] ⋅ uband[i].vectors[:, bandindex]
# 			W *= uu
# 		end
# 		uu = sum(i -> cis(2π * ([1, 0, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, bandindex]) * uband[end].vectors[i, bandindex], eachindex(TB.orb_location))
# 		W = uu * W
# 		W = W / abs(W)
# 	elseif bandindex isa AbstractVector{<:Integer}
# 		W = I
# 		for i in 1:Nk-1
# 			uu = [uband[i+1].vectors[:, m] ⋅ uband[i].vectors[:, n] for m in bandindex, n in bandindex]
# 			W = uu * W
# 		end
# 	else
# 		error("Wrong bandindex.")
# 	end

# end

