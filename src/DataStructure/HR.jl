export HR, numorb, shift_fermi_energy!, spinHR, reindexHR, reindexHR!
"""
The data in wannier90_hr.dat.
"""
struct HR{T <: Number}
	orbindex::Vector{Int}
	path::Matrix{Int} # [N,5]
	value::Vector{T}
end
function HR(path::AbstractMatrix{<:Integer}, value::AbstractVector{T};
	orbindex = sort(unique(path[:, 4:5])), μ::Real = 0, hrsort = 'N')::HR{T} where {T <: Number}

	hr = HR{T}(orbindex, path, value)

	if abs(μ) > 1e-6
		shift_fermi_energy!(hr, μ)
	end

	if hrsort[1] ∈ ['Y', 'y']
		sort!(hr)
	end

	return hr
end

function Base.show(io::IO, hr::HR)
	print(io, "HR with $(length(hr)) hoppings and $(numorb(hr)) orbitals.")
end
Base.length(hr::HR) = length(hr.value)
@inline numorb(hr::HR) = length(hr.orbindex)

Base.convert(::HR{T₁}, hr::HR{T₂}) where {T₁, T₂} =
	HR{T₁}(hr.orbindex, hr.path, convert(Vector{T₁}, hr.value))
Base.convert(::Type{HR{T₁}}, hr::HR{T₂}) where {T₁ <: Number, T₂} =
	HR{T₁}(hr.orbindex, hr.path, convert(Vector{T₁}, hr.value))

"""
	filter(F, hr::HR{T}) where {T}

F(i, j, R, value) -> Bool.
t_{i0,jR} = value.
"""
function Base.filter(F, hr::HR{T}) where {T}

	I = similar(hr.value, Bool)
	for i in eachindex(I)
		I[i] = F(hr.path[i, 4], hr.path[i, 5], hr.path[i, 1:3], hr.value[i])
	end

	new_path = hr.path[I, :]
	new_value = hr.value[I]

	return HR{T}(hr.orbindex, new_path, new_value)
end

"""
	sort!(hr::HR{T}) where {T}

Sort hr.path and hr.value by hr.path.
"""
function Base.sort!(hr::HR{T}) where {T}
	pathI = sortslices([hr.path collect(eachindex(hr.value))]; dims = 1)
	I = pathI[:, end]

	hr.path .= pathI[:, 1:5]

	value_copy = copy(hr.value)
	hr.value .= value_copy[I]

	return hr
end

function shift_fermi_energy!(hr::HR, μ::Real)
	I = findall(x -> iszero(x[1:3]) && isequal(x[4], x[5]), eachrow(hr.path))
	hr.value[I] .-= μ
	return hr
end

"""
	spinHR(hr::HR{T}; buildindex = 'Y')::HR{T} where {T}

Generate a spinful hr from a spinless hr.
Acctually, `spinHR` just makes a copy, it won't introduce new hopping.
"""
function spinHR(hr::HR{T}) where {T}
	up_path = hr.path
	up_value = hr.value

	norb = length(hr.orbindex)
	#make sure using a seperate path
	dn_path = copy(hr.path)
	dn_path[:, 4:5] .+= norb
	#conj is fron time-reversal.
	dn_value = conj.(hr.value)

	sum_path = [up_path; dn_path]
	sum_value = [up_value; dn_value]

	sum_orbindex = sort([hr.orbindex; hr.orbindex .+ norb])

	return HR{T}(sum_orbindex, sum_path, sum_value)
end



















"""
	reindexHR!(hr::HR, newindex::AbstractVector{<:Integer})::HR

This function will change the orbindex and orbpath of hr, and don't need the index of hr.
"""
function reindexHR(hr::HR{T}, newindex::AbstractVector{<:Integer}; buildindex = 'Y')::HR{T} where {T}

	if numorb(hr) ≠ length(newindex)
		error("Wrong length of newindex to hr.")
	end

	neworbpath = Matrix{Int}(undef, length(hr), 2)
	@views for i in eachindex(newindex)
		I = hr.path[:, 4:5] .== hr.orbindex[i]
		neworbpath[I] .= newindex[i]
	end

	newpath = [hr.path[:, 1:3] neworbpath]

	return HR{T}(newpath, deepcopy(hr.value); buildindex)
end
"""
	reindexHR(hr::HR, N::Integer)::HR
	reindexHR!(hr::HR, N::Integer)::HR

This function will change the orbindex and orbpath of hr, and add N to these two data.
"""
function reindexHR(hr::HR{T}, N::Integer)::HR{T} where {T}
	return reindexHR!(HR{T}(deepcopy(hr.path), deepcopy(hr.value); buildindex = 'N'), N)
end
function reindexHR!(hr::HR{T}, N::Integer)::HR{T} where {T}

	hr.path[:, 4:5] .+= N

	return hr
end
