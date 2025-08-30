export unstackPOSCAR
"""
	unstackPOSCAR(poscar::POSCAR, zlist::Real...; unit = "Cartesian", zindex::Integer = 3)

Only used to split a multilayer POSCAR to obtain each layer's POSCAR.
`poscar` is the multilayer POSCAR;
`zlist` is the list of interface between adjacent layers;
`unit` is the unit of `zlist`, can set as "Cartesian" or "Direct";
`zindex` represents the direction of out-of-plane, can set as 1, 2 or 3.
"""
function unstackPOSCAR(poscar::POSCAR, zlist::Real...; unit = "Cartesian", zindex::Integer = 3)

	if unit[1] ∈ ['C', 'c']
	elseif unit[1] ∈ ['D', 'd']
		zlist = zlist .* poscar.lattvec[zindex, zindex]
	else
		error("Wrong unit from unstackPOSCAR.")
	end

	atomzlocation = [x[zindex] for x in poscar.location]

	zlist = [minimum(atomzlocation) - 1; sort([zlist...])]

	n = length(zlist)
	layerposcar = Vector{POSCAR}(undef, n)

	for i in 1:n-1
		I = findall(z -> zlist[i] < z < zlist[i+1], atomzlocation)
		layerposcar[i] = POSCAR(
			poscar.lattvec,
			poscar.location[I];
			name = poscar.name[I],
			index = poscar.index[I],
			periodicity = poscar.periodicity,
		)
	end
	I = atomzlocation .> zlist[n]
	layerposcar[n] = POSCAR(
		poscar.lattvec,
		poscar.location[I];
		name = poscar.name[I],
		index = poscar.index[I],
		periodicity = poscar.periodicity,
	)

	if poscar.num ≠ sum(p -> p.num, layerposcar)
		error("Some atoms might happen to be at the interface.")
	end

	return layerposcar
end

function findinterface(poscar::POSCAR; zindex::Integer = 3, Δ::Real = 3)

	distribute = GenerateDistribute(poscar, zindex; σ = Δ / 8)

	distribute1 = diff(distribute)
	distribute2 = diff(distribute1)

	I1 = abs.(distribute1) .< 0.1

	return zlist
end

function AMPD(data::AbstractArray{T}) where {T <: Real}
    """
    实现AMPD算法
    :param data: 1-D AbstractArray{T} 
    :return: 波峰所在索引值的数组
    """
    p_data = zeros(Int, length(data))
    count = length(data)
    arr_rowsum = Int[]
    for k in 1:div(count, 2)
        row_sum = 0
        for i in k:(count - k)
            if data[i] > data[i - k] && data[i] > data[i + k]
                row_sum -= 1
            end
        end
        push!(arr_rowsum, row_sum)
    end
    min_index = findmin(arr_rowsum)[2]
    max_window_length = min_index
    for k in 1:max_window_length
        for i in k:(count - k)
            if data[i] > data[i - k] && data[i] > data[i + k]
                p_data[i] += 1
            end
        end
    end
    return findall(x -> x == max_window_length, p_data)
end

function GenerateDistribute(poscar::POSCAR, axis::Integer; σ::Real = 0.5)

	axislocation = map(r -> r[axis], poscar.location)
	axis = poscar.lattvec[axis, axis]
	axis = range(-axis / 2, 3 * axis / 2, step = 0.05)

	axisdensity = zeros(length(axis))
	for b in axislocation
		GF = Gaussian1D(; a = 1, b, c = σ)
		axisdensity .+= GF.(axis)
	end

	return axisdensity
end
