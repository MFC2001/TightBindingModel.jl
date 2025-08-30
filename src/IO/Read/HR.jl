export ReadHR

"""
	ReadHR(HRfile::AbstractString, heps::Real=0; fermienergy::Real=0)::HR

Read wannier90_hr.dat, return HR.
"""
function ReadHR(file::AbstractString, heps::Real = 0; readimag = 'Y', μ::Real = 0, hrsort = 'N')::HR

	if heps < 0
		error("Please input a positive heps!")
	end

	file = open(file, "r")
	readline(file)
	norb = parse(Int, readline(file))

	for _ in 1:((parse(Int, readline(file))-1)÷15+1)
		readline(file)
	end
	hr = readdlm(file)
	close(file)


	allpath = Int.(hr[:, 1:5])
	allvalue = hr[:, 6:7]

	if readimag[1] ∈ ['N', 'n'] || iszero(allvalue[:, 2])
		allvalue = allvalue[:, 1]
		I = map(x -> abs(x) >= heps, allvalue)

		path = allpath[I, :]
		value = allvalue[I]
	else
		allvalue = complex.(allvalue[:, 1], allvalue[:, 2])
		heps2 = heps^2
		I = map(x -> abs2(x) >= heps2, allvalue)

		path = allpath[I, :]
		value = allvalue[I]
	end


	return HR(path, value; μ, hrsort)
end
