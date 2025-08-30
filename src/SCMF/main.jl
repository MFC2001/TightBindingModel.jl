export SCMF
"""

"""
function SCMF(
	kgrid::MonkhorstPack,
	hr::HR,
	poscar::POSCAR,
	orbital::ORBITAL;
	U::Function,
	J::Function,
	Ne = Int(length(orbital) / 2),
	T = 0,
	Δ₀ = nothing,
	Δeps = 1e-8,
	heps = 1e-6,
)

	irredkgrid = irreducible_kcoords(kgrid, poscar; symtol = 1e-5, check_symmetry = true)
	get_lattsymop!(irredkgrid, poscar, orbital)


	H₀ = Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}(undef, length(irredkgrid))
	Threads.@threads for i in eachindex(irredkgrid.kcoord)
		H₀[i] = Hamilton(irredkgrid.kcoord[i], hr, poscar.basis)
	end


	norb = length(orbital)
	if T == 0
		T = zeros(Bool, norb)
		T[1:Ne] .= 1
		FD(x) = T
	else
		error("To be continued for T ≠ 0.")
	end
	Vik = coulomb_inner_product(irredkgrid, U, J)
	DM = DensityMatrix(irredkgrid, FD)

	ijkittr = [(i, j, k) for k in eachindex(irredkgrid.redkcoord), j in 1:norb, i in 1:norb]
	order(band) = MForder!(Δ, Δittr, ijkittr, Vik, DM(band))

	if isnothing(Δ₀)
		# band = BAND(Hk)
		# Δ₀ = order(Vk, DensityMatrix(FD, band), minusmap)

		Δ₀ = [rand(ComplexF64, norb, norb) for _ in 1:length(irredkgrid)]
		scale = maximum(abs, H₀[1]) / maximum(x -> maximum(abs, x), Δ₀)
		Δ₀ = Δ₀ .* scale
	end


	Hk = SCMF_kernal(H₀, order, Δ₀, Δeps)

	return HK2HR(Hk, irredkgrid, poscar.basis, heps)
end

