
#本质上，要求输入的q能忠实地反映kgrid到自身的映射。
#要求任意的k+q属于kgrid，自然k_Γ+q也属于k_Γ。
#所有满足条件的q均隶属于一个k+k组成的集合，这个k点集合恰好就是kgrid_Γ。
#当输入的q隶属于一个网格间距为kgrid整数倍的qgrid时，若倍数为奇数则要求无偏移，若倍数为偶数则可以有偏移，此时q也一定隶属于kgrid_Γ.
#事实上，目前该数据结构实现的优化仅有：避免了bandkq与Jkq的重复计算。

struct BSEqgrid{
	HT <: Number,
	HU <: ReciprocalHoppings{HT},
	TBT <: TightBindModel{HT, HU},
	U_d <: Number,
	LU_d <: GaussLRCorrection,
	J¹_d <: Number,
	J²_d <: Number,
	U_x <: Number,
	LU_x <: GaussLRCorrection,
	J¹_x <: Number,
	J²_x <: Number,
} <: AbstractBSE
	TB::TBT
	type::Symbol
	scissor::Float64
	kgrid::RedKgrid
	vckmap::vckMap
	kgrid_Γ::RedKgrid
	#kgrid to kgrid_Γ
	addmap::Matrix{Int}
	minusmap::Matrix{Int}
	#for kgrid.
	bandk::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	#for kgrid_Γ.
	Kᵈ_U::UwithLR{U_d, LU_d} #目前UwithLR未包含并行优化
	Kᵈ_J¹::BaseReciprocalHoppings{J¹_d}
	Kᵈ_J²::BaseReciprocalHoppings{J²_d}
	Kˣ_U::UwithLR{U_x, LU_x}
	Kˣ_J¹::BaseReciprocalHoppings{J¹_x}
	Kˣ_J²::BaseReciprocalHoppings{J²_x}
	#pre_calculated
	Kᵈ_Uk::Vector{Matrix{ComplexF64}}
	Kᵈ_J²k::Vector{Matrix{ComplexF64}}
	Kˣ_J¹k::Vector{Matrix{ComplexF64}}
	Kˣ_J²k::Vector{Matrix{ComplexF64}}
	#preallocated buffer
	Kᵈ_J¹q::Matrix{ComplexF64}
	Kˣ_Uq::Matrix{ComplexF64}
end
function Base.show(io::IO, bse::BSEqgrid)
	print(io, "$(ndims(bse.TB)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end

"""
"""
function BSEqgrid(TB::AbstractTightBindModel, kgrid::MonkhorstPack, v, c, type::Symbol;
	scissor::Real = 0,
	Kᵈ_U_LR::Bool = true,
	Kˣ_U_LR::Bool = true,
	αrcut::Real = 4.5,
	rcut::Real = 10,
	ϵ = 1,
	δ::Real = 1e-6,
	Kᵈ_U::Union{HR, <:AbstractString, Nothing} = nothing,
	Kᵈ_J¹::Union{HR, <:AbstractString, Nothing} = nothing,
	Kᵈ_J²::Union{HR, <:AbstractString, Nothing} = nothing,
	Kˣ_U::Union{HR, <:AbstractString, Nothing} = nothing,
	Kˣ_J¹::Union{HR, <:AbstractString, Nothing} = nothing,
	Kˣ_J²::Union{HR, <:AbstractString, Nothing} = nothing,
)
	if type ∉ [:NP, :SP]
		error("BSE type should be :NP or :SP!")
	end

	kgrid = RedKgrid(kgrid)
	vckmap = vckMap(v, c, length(kgrid))
	bandk = BAND(kgrid, TB; vector = true)
	_sum_wave_is_real!.(bandk)

	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	addmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, +)
	minusmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, -)


	norb = numorb(TB)
	Kᵈ_J¹ = _BSE_preprocess_J(Kᵈ_J¹, norb)
	Kᵈ_J² = _BSE_preprocess_J(Kᵈ_J², norb)
	Kˣ_J¹ = _BSE_preprocess_J(Kˣ_J¹, norb)
	Kˣ_J² = _BSE_preprocess_J(Kˣ_J², norb)


	α = αrcut / rcut
	Kᵈ_U = _BSE_preprocess_U(Kᵈ_U, Kᵈ_U_LR, TB, α, ϵ, δ)
	Kˣ_U = _BSE_preprocess_U(Kˣ_U, Kˣ_U_LR, TB, α, 1, δ)


	#preprocessing
	nk = length(kgrid)
	Kᵈ_Uk = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Kᵈ_J²k = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Kˣ_J¹k = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Kˣ_J²k = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	Threads.@threads for k in 1:nk
		Kᵈ_U(Kᵈ_Uk[k], kgrid_Γ[k], nk)
		Kᵈ_Uk[k] ./= nk
		Kᵈ_J²(Kᵈ_J²k[k], kgrid_Γ[k], nk)
		Kᵈ_J²k[k] ./= nk
		Kˣ_J¹(Kˣ_J¹k[k], kgrid_Γ[k])
		Kˣ_J¹k[k] ./= nk
		Kˣ_J²(Kˣ_J²k[k], kgrid_Γ[k])
		Kˣ_J²k[k] ./= nk
	end

	#preallocated buffer
	Kᵈ_J¹q = Matrix{ComplexF64}(undef, norb, norb)
	Kˣ_Uq = Matrix{ComplexF64}(undef, norb, norb)
	return BSEqgrid(TB, type, Float64(scissor), kgrid, vckmap, kgrid_Γ, addmap, minusmap, bandk,
		Kᵈ_U, Kᵈ_J¹, Kᵈ_J², Kˣ_U, Kˣ_J¹, Kˣ_J²,
		Kᵈ_Uk, Kᵈ_J²k, Kˣ_J¹k, Kˣ_J²k, Kᵈ_J¹q, Kˣ_Uq)
end

