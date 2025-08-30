export HK2HR
"""


"""
function HK2HR(Hk::AbstractArray{<:Number}, kgrid::RedKgrid, heps::Real = 1e-6)

	Nk = length(kgrid)
	Rgrid = map(k -> Vec3{Int}(k .* kgrid.kgrid_size), kgrid.kdirect)

	norb = size(Hk, 1)

	heps2 = heps^2

	path = Matrix{Int}(undef, 5, 0)
	value = Vector{ComplexF64}(undef, 0)

	lk = ReentrantLock()
	Threads.@threads for aimpath in Rgrid

		t = zeros(ComplexF64, norb, norb)
		for k in 1:Nk
			t .+= Hk[:, :, k] .* cis(-2π * (kgrid.kdirect[k] ⋅ aimpath))
		end
		t ./= Nk

		I = abs2.(t) .> heps2
		value_t = t[I]

		II = Tuple.(CartesianIndices(t)[I])
		nt = length(II)
		path_t = Matrix{Int}(undef, 5, nt)
		for (i, IIi) in enumerate(II)
			path_t[:, i] = [aimpath; IIi[1]; IIi[2]]
		end

		lock(lk) do
			path = [path path_t]
			append!(value, value_t)
		end
	end

	return HR(path', value; hrsort = "Y")
end

function HK2HR(Hk_irred::AbstractVector{<:AbstractMatrix{<:Number}}, irredkgrid::IrredKgrid, heps::Real = 1e-6)

	Nk = length(irredkgrid.redkcoord)
	norb = size(Hk_irred[1], 1)
	gridindex = map(k -> Vec3{Int}(k .* irredkgrid.kgrid_size), irredkgrid.redkdirect)

	Hk = Vector{eltype(Hk_irred)}(undef, Nk)
	Threads.@threads for k in eachindex(Hk)
		U = irredkgrid.lattsymop[k].U
		Hk[k] = Hermitian(U * Hk_irred[irredkgrid.irmap[k]] * U')
	end


	heps2 = heps^2

	path = Matrix{Int}(undef, 0, 5)
	value = Vector{ComplexF64}(undef, 0)

	lk = ReentrantLock()
	Threads.@threads for aimpath in gridindex

		t = zeros(ComplexF64, norb, norb)
		for k in 1:Nk
			t .+= Hk[k][:, :] .* cis(-2π * (irredkgrid.redkcoord[k] ⋅ aimpath))
		end
		t ./= Nk

		I = abs2.(t) .> heps2
		value_t = t[I]

		I = Tuple.(CartesianIndices(t)[I])
		nt = length(I)
		path_t = Matrix{Int}(undef, 5, nt)
		for (i, I) in enumerate(I)
			path_t[:, i] = [aimpath..., I[1], I[2]]
		end

		lock(lk) do
			path = [path; transpose(path_t)]
			append!(value, value_t)
		end
	end

	return HR(path, value; hrsort = "Y")
end

# function getucpathfromkmesh(counter, kmesh)

# 	meshindex = kmesh.meshindex

# 	T = Set([-counter, counter])
# 	I = map(x -> x[1] ∈ T || x[2] ∈ T || x[3] ∈ T, meshindex)

# 	meshindex = meshindex[I, :]

# 	I = map(x -> abs(x[1]) <= counter && abs(x[2]) <= counter && abs(x[3]) <= counter, meshindex)

# 	return meshindex[I, :]
# end
