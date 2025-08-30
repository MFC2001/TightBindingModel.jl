
function Kernal(bse::BSE, q::ReducedCoordinates)
	# We can't calculate Γ directly.
	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)
	end

	bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
	_sum_wave_is_real!.(bandkq)

	nk = length(bse.kgrid)
	norb = numorb(bse.TB)

	#the cost of calculating J_matrix is small.
	Kᵈ_J²kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Kˣ_J²kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Threads.@threads for k in 1:nk
		kq = bse.kgrid_Γ[k] + q
		bse.Kᵈ_J²(Kᵈ_J²kq[k], kq)
		bse.Kˣ_J²(Kˣ_J²kq[k], kq)
		Kᵈ_J²kq[k] ./= nk
		Kˣ_J²kq[k] ./= nk
	end

	Kᵈ_J¹q = bse.Kᵈ_J¹(q) ./ nk
	#May parallel optimization is required here.
	Kˣ_Uq = bse.Kˣ_U(q) ./ nk

	Kᵈ = CreatKᵈ_Q(bse.bandk, bandkq, bse.Kᵈ_Uk, Kᵈ_J¹q, Kᵈ_J²kq, bse.addmap, bse.minusmap)
	Kˣ = CreatKˣ_Q(bse.bandk, bandkq, Kˣ_Uq, bse.Kˣ_J¹k, Kˣ_J²kq, bse.addmap, bse.minusmap)

	return bandkq, Kᵈ, Kˣ
end
function Kernal!(bse::BSE, q::AbstractVector)
	#Make sure you know that you can't calculate more than one Hbse at the same time.

	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)
	end

	bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
	_sum_wave_is_real!.(bandkq)

	nk = length(bse.kgrid)

	#the cost of calculating J_matrix is small.
	Threads.@threads for k in 1:nk
		kq = bse.kgrid_Γ[k] + q
		bse.Kᵈ_J²(bse.Kᵈ_J²kq[k], kq)
		bse.Kˣ_J²(bse.Kˣ_J²kq[k], kq)
		bse.Kᵈ_J²kq[k] ./= nk
		bse.Kˣ_J²kq[k] ./= nk
	end

	bse.Kᵈ_J¹(bse.Kᵈ_J¹q, q)
	bse.Kᵈ_J¹q ./= nk
	#May parallel optimization is required here.
	bse.Kˣ_U(bse.Kˣ_Uq, q)
	bse.Kˣ_Uq ./= nk

	Kᵈ = CreatKᵈ_Q(bse.bandk, bandkq, bse.Kᵈ_Uk, bse.Kᵈ_J¹q, bse.Kᵈ_J²kq, bse.addmap, bse.minusmap)
	Kˣ = CreatKˣ_Q(bse.bandk, bandkq, bse.Kˣ_Uq, bse.Kˣ_J¹k, bse.Kˣ_J²kq, bse.addmap, bse.minusmap)

	return bandkq, Kᵈ, Kˣ
end

function Kernal(bse::BSEqgrid, q::AbstractVector)

	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)

		bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
		_sum_wave_is_real!.(bandkq)

		nk = length(bse.kgrid)
		norb = numorb(bse.TB)

		#the cost of calculating J_matrix is small.
		Kᵈ_J²kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
		Kˣ_J²kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
		Threads.@threads for k in 1:nk
			kq = bse.kgrid_Γ[k] + q
			bse.Kᵈ_J²(Kᵈ_J²kq[k], kq)
			bse.Kˣ_J²(Kˣ_J²kq[k], kq)
			Kᵈ_J²kq[k] ./= nk
			Kˣ_J²kq[k] ./= nk
		end

		Kᵈ_J¹q = bse.Kᵈ_J¹(q) ./ nk
		#May parallel optimization is required here.
		Kˣ_Uq = bse.Kˣ_U(q) ./ nk

	else
		nk = length(bse.kgrid)

		kq_kindex = Vector{Int}(undef, nk)
		kΓq_kindex = Vector{Int}(undef, nk)
		Threads.@threads for ki in 1:nk
			kq = bse.kgrid[ki] + q
			kq_kindex[ki] = findfirst(k -> all(isinteger, k - kq), bse.kgrid)
			kΓq = bse.kgrid_Γ[ki] + q
			kΓq_kindex[ki] = findfirst(k -> all(isinteger, k - kΓq), bse.kgrid_Γ)
		end

		bandkq = bse.bandk[kq_kindex]
		Kᵈ_J²kq = bse.Kᵈ_J²k[kΓq_kindex]
		Kˣ_J²kq = bse.Kˣ_J²k[kΓq_kindex]

		Kᵈ_J¹q = bse.Kᵈ_J¹(q) ./ nk
		#May parallel optimization is required here.
		Kˣ_Uq = bse.Kˣ_U(q) ./ nk
	end

	Kᵈ = CreatKᵈ_Q(bse.bandk, bandkq, bse.Kᵈ_Uk, Kᵈ_J¹q, Kᵈ_J²kq, bse.addmap, bse.minusmap)
	Kˣ = CreatKˣ_Q(bse.bandk, bandkq, Kˣ_Uq, bse.Kˣ_J¹k, Kˣ_J²kq, bse.addmap, bse.minusmap)

	return bandkq, Kᵈ, Kˣ
end
function Kernal!(bse::BSEqgrid, q::AbstractVector)

	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)

		bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
		_sum_wave_is_real!.(bandkq)

		nk = length(bse.kgrid)
		norb = numorb(bse.TB)

		#the cost of calculating J_matrix is small.
		Kᵈ_J²kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
		Kˣ_J²kq = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
		Threads.@threads for k in 1:nk
			kq = bse.kgrid_Γ[k] + q
			bse.Kᵈ_J²(Kᵈ_J²kq[k], kq)
			bse.Kˣ_J²(Kˣ_J²kq[k], kq)
			Kᵈ_J²kq[k] ./= nk
			Kˣ_J²kq[k] ./= nk
		end

		bse.Kᵈ_J¹(bse.Kᵈ_J¹q, q)
		bse.Kᵈ_J¹q ./= nk
		#May parallel optimization is required here.
		bse.Kˣ_U(bse.Kˣ_Uq, q)
		bse.Kˣ_Uq ./= nk

	else
		nk = length(bse.kgrid)

		kq_kindex = Vector{Int}(undef, nk)
		kΓq_kindex = Vector{Int}(undef, nk)
		Threads.@threads for ki in 1:nk
			kq = bse.kgrid[ki] + q
			kq_kindex[ki] = findfirst(k -> all(isinteger, k - kq), bse.kgrid)
			kΓq = bse.kgrid_Γ[ki] + q
			kΓq_kindex[ki] = findfirst(k -> all(isinteger, k - kΓq), bse.kgrid_Γ)
		end

		bandkq = bse.bandk[kq_kindex]
		Kᵈ_J²kq = bse.Kᵈ_J²k[kΓq_kindex]
		Kˣ_J²kq = bse.Kˣ_J²k[kΓq_kindex]

		bse.Kᵈ_J¹(bse.Kᵈ_J¹q, q)
		bse.Kᵈ_J¹q ./= nk
		#May parallel optimization is required here.
		bse.Kˣ_U(bse.Kˣ_Uq, q)
		bse.Kˣ_Uq ./= nk
	end

	Kᵈ = CreatKᵈ_Q(bse.bandk, bandkq, bse.Kᵈ_Uk, bse.Kᵈ_J¹q, bse.Kᵈ_J²kq, bse.addmap, bse.minusmap)
	Kˣ = CreatKˣ_Q(bse.bandk, bandkq, bse.Kˣ_Uq, bse.Kˣ_J¹k, bse.Kˣ_J²kq, bse.addmap, bse.minusmap)

	return bandkq, Kᵈ, Kˣ
end
function CreatKᵈ_Q(bandk, bandkq, Uk, Jq, Jkq, addmap, minusmap)

	norb = size(bandk[1].vectors, 1)
	index = Tuple.(CartesianIndices((norb, norb)))

	kernal = function (v′, c′, k′, v, c, k)
		ψ₁ = view(bandkq[k′].vectors, :, c′)
		ψ₂ = view(bandk[k].vectors, :, v)
		ψ₃ = view(bandk[k′].vectors, :, v′)
		ψ₄ = view(bandkq[k].vectors, :, c)
		T₁ = conj.(ψ₁ * transpose(ψ₂))
		T₂ = ψ₃ * transpose(ψ₄)
		U₁ = Uk[minusmap[k′, k]]
		J₃ = Jkq[addmap[k′, k]]
		# note this minus
		return -sum(index) do (i, j)
			T₁[i, j] * T₂[j, i] * U₁[i, j] + T₁[i, j] * T₂[i, j] * Jq[i, j] + T₁[i, i] * T₂[j, j] * J₃[i, j]
		end
	end

	return kernal
end
function CreatKˣ_Q(bandk, bandkq, Uq, Jk, Jkq, addmap, minusmap)

	norb = size(bandk[1].vectors, 1)
	index = Tuple.(CartesianIndices((norb, norb)))

	kernal = function (v′, c′, k′, v, c, k)
		ψ₁ = view(bandkq[k′].vectors, :, c′)
		ψ₂ = view(bandk[k].vectors, :, v)
		ψ₃ = view(bandkq[k].vectors, :, c)
		ψ₄ = view(bandk[k′].vectors, :, v′)
		T₁ = conj.(ψ₁ * transpose(ψ₂))
		T₂ = ψ₃ * transpose(ψ₄)
		J₂ = Jk[minusmap[k′, k]]
		J₃ = Jkq[addmap[k′, k]]
		return sum(index) do (i, j)
			T₁[i, j] * T₂[j, i] * Uq[i, j] + T₁[i, j] * T₂[i, j] * J₂[i, j] + T₁[i, i] * T₂[j, j] * J₃[i, j]
		end
	end

	return kernal
end
