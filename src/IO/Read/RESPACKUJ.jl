export ReadRESPACKUJ
"""
	Read dat.W(JVX)mat or dat.h_mat_r from RESPACK.
"""
function ReadRESPACKUJ(file::AbstractString, norb::Integer; onsite = "Y")

	norb2 = norb^2

	file = open(file, "r")
	for _ in 1:3
		readline(file)
	end

	ppath = Vector{Vector{Int}}(undef, 0)
	vvalue = Vector{Matrix{Float64}}(undef, 0)
	while !eof(file)
		mark(file)
		line = readline(file)
		if occursin(r"^\s*(#|$)", line)
			continue
		end

		push!(ppath, parse.(Int, split(line)))

		M = Matrix{Float64}(undef, norb2, 4)
		for i in 1:norb2
			M[i, :] = parse.(Float64, split(readline(file)))
		end
		push!(vvalue, M)
	end

	close(file)

	n = length(ppath)
	path = Matrix{Int}(undef, n * norb2, 5)
	value = Vector{ComplexF64}(undef, n * norb2)
	@views if onsite[1] ∈ ['Y', 'y']
		for i in 1:n, (ii, j) in enumerate(norb2*(i-1)+1:norb2*i)
			path[j, :] = [ppath[i]; vvalue[i][ii, 1:2]]
			value[j] = complex.(vvalue[i][ii, 3], vvalue[i][ii, 4])
		end
	elseif onsite[1] ∈ ['N', 'n']
		for i in 1:n
			if iszero(ppath[i])
				for (ii, j) in enumerate(norb2*(i-1)+1:norb2*i)
					path[j, :] .= [ppath[i]; vvalue[i][ii, 1:2]]
					if path[j, 4] == path[j, 5]
						value[j] = 0
					else
						value[j] = complex.(vvalue[i][ii, 3], vvalue[i][ii, 4])
					end
				end
			else
				for (ii, j) in enumerate(norb2*(i-1)+1:norb2*i)
					path[j, :] = [ppath[i]; vvalue[i][ii, 1:2]]
					value[j] = complex.(vvalue[i][ii, 3], vvalue[i][ii, 4])
				end
			end
		end
	else
		error("Wrong keyword onsite.")
	end


	return HR(path, value; hrsort = 'N', buildhop = 'Y')
end
