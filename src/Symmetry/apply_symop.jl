"""
Apply a symmetry operation to eigenvectors `ψk` at a given `kpoint` to obtain an
equivalent point in [-0.5, 0.5)^3 and associated eigenvectors (expressed in the
basis of the new ``k``-point).
"""
#symop::LattSymOp
function apply_symop(symop::SymOp, basis, kpoint, ψk::AbstractVecOrMat)
	S, τ = symop.S, symop.τ
	isone(symop) && return kpoint, ψk

	# Apply S and reduce coordinates to interval [-0.5, 0.5)
	# Doing this reduction is important because
	# of the finite kinetic energy basis cutoff
	@assert all(-0.5 .≤ kpoint.coordinate .< 0.5)
	Sk_raw = S * kpoint.coordinate
	Sk = normalize_kpoint_coordinate(Sk_raw)
	kshift = convert.(Int, Sk - Sk_raw)
	@assert all(-0.5 .≤ Sk .< 0.5)

	# Check whether the resulting k-point is in the basis:
	ikfull = findfirst(1:length(basis.kpoints)) do idx
		all(isinteger, basis.kpoints[idx].coordinate - Sk)
	end
	if isnothing(ikfull)
		# Build new k-point datastructure
		Skpoint = Kpoint(basis, Sk, kpoint.spin)
	else
		Skpoint = basis.kpoints[ikfull]
		@assert Skpoint.coordinate ≈ Sk
	end

	# Since the eigenfunctions of the Hamiltonian at k and Sk satisfy
	#      u_{Sk}(x) = u_{k}(S^{-1} (x - τ))
	# their Fourier transform is related as
	#      ̂u_{Sk}(G) = e^{-i G \cdot τ} ̂u_k(S^{-1} G)
	invS = Mat3{Int}(inv(S))
	Gs_full = [G + kshift for G in G_vectors(basis, Skpoint)]
	ψSk = zero(ψk)
	for iband ∈ axes(ψk, 2)
		for (ig, G_full) in enumerate(Gs_full)
			igired = index_G_vectors(basis, kpoint, invS * G_full)
			@assert igired !== nothing
			ψSk[ig, iband] = cis2pi(-dot(G_full, τ)) * ψk[igired, iband]
		end
	end

	Skpoint, ψSk
end