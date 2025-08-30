export ContinuationHR
"""

F(i,j,R) = tᵢⱼ(R), e.g. i0 and jR, R is reduced coordinate.
Pay attention to the bonds that exceed periodic boundary.
"""

function ContinuationHR(hr::HR, F::Function, cell::Cell, maxlattdist::Real)
	ucpath = gridindex(cell, maxlattdist)
	return ContinuationHR(hr, F, ucpath)
end
function ContinuationHR(hr::HR, F::Function, grid_size::AbstractVector{<:Integer})
	ucpath = gridindex(grid_size)
	return ContinuationHR(hr, F, ucpath)
end
function ContinuationHR(hr::HR, F::Function, ucpath::AbstractVector{<:AbstractVector})

	hr_path = Matrix{Int}(undef, 5, 0)
	hr_value = Vector{typeof(F(1, 1, [1, 1, 1]))}(undef, 0)

	lk = ReentrantLock()
	norb = numorb(hr)
	Threads.@threads for I in CartesianIndices((norb, norb))

		existpath = path.(hr.hop[I])
		addpathindex = findall(x -> x ∉ existpath, ucpath)

		(i, j) = Tuple(I)
		addvalue = map(r -> F(i, j, ucpath[r]), addpathindex)


		addpath = Matrix{Int}(undef, 5, length(addvalue))
		for ii in eachindex(addvalue)
			addpath[:, ii] = [ucpath[addpathindex[ii]]..., i, j]
		end

		lock(lk) do
			hr_path = [hr_path addpath]
			append!(hr_value, addvalue)
		end
	end

	hr_path = [hr.path; transpose(hr_path)]
	hr_value = [hr.value; hr_value]

	return HR(hr_path, hr_value; hrsort = 'Y')
end
function ContinuationHR(hr::HR, F::Function, ucpath::AbstractMatrix{<:AbstractVector{<:AbstractVector}})

	hr_path = Matrix{Int}(undef, 5, 0)
	hr_value = Vector{typeof(F(1, 1, [1, 1, 1]))}(undef, 0)

	lk = ReentrantLock()
	norb = numorb(hr)
	Threads.@threads for I in CartesianIndices((norb, norb))

		existpath = path.(hr.hop[I])
		addpathindex = findall(x -> x ∉ existpath, ucpath[I])

		(i, j) = Tuple(I)
		addvalue = map(r -> F(i, j, ucpath[I][r]), addpathindex)


		addpath = Matrix{Int}(undef, 5, length(addvalue))
		for ii in eachindex(addvalue)
			addpath[:, ii] = [ucpath[I][addpathindex[ii]]..., i, j]
		end

		lock(lk) do
			hr_path = [hr_path addpath]
			append!(hr_value, addvalue)
		end
	end

	hr_path = [hr.path; transpose(hr_path)]
	hr_value = [hr.value; hr_value]

	return HR(hr_path, hr_value; hrsort = 'Y')
end
