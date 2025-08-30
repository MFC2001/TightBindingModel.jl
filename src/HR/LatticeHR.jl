export LatticeHR
"""

F(i,j,R) = tᵢⱼ(R), e.g. i0 and jR, R is reduced coordinate.
Pay attention to the bonds that exceed periodic boundary.
"""
function LatticeHR(F, grid::MonkhorstPack, norb::Integer)
	ucpath = gridindex(grid.kgrid_size)
	return LatticeHR(F, ucpath, norb)
end
function LatticeHR(F, grid::AbstractVector{<:Integer}, norb::Integer)
	ucpath = gridindex(grid)
	return LatticeHR(F, ucpath, norb)
end
function LatticeHR(F, ucpath::AbstractVector{<:AbstractVector}, norb::Integer)

	hr_path = Matrix{Int}(undef, 5, 0)
	hr_value = Vector{Float64}(undef, 0)

	lk = ReentrantLock()
	Threads.@threads for I in CartesianIndices((norb, norb))

		(i, j) = Tuple(I)
		addvalue = map(R -> F(i, j, R), ucpath)


		addpath = Matrix{Int}(undef, 5, length(addvalue))
		for ii in eachindex(addvalue)
			addpath[:, ii] = [ucpath[ii]..., i, j]
		end

		lock(lk) do
			hr_path = [hr_path addpath]
			hr_value = [hr_value; addvalue]
		end
	end


	return HR(transpose(hr_path), hr_value; orbindex = collect(1:norb), hrsort = 'Y')
end
