export integrate
"""
	integrateHR(hrs::HR...; hrsort = "Y", buildindex = "Y")::HR

P.S. This function will try to create a new hr with the same type of hrs[1], so make sure that all hrs satisfy the condition.
"""
function integrate(hrs::HR...; hrsort = "Y", buildindex = "Y")::HR

	path = Vector{eltype(hrs[1].path)}(undef, 0)
	value = Vector{eltype(hrs[1].value)}(undef, 0)
	orbpath = Matrix{eltype(hrs[1].orbpath)}(undef, 2, 0)
	for hr in hrs
		append!(path, hr.path)
		append!(value, hr.value)
		orbpath = [orbpath hr.orbpath]
	end

	return HR(path, orbpath, value; hrsort, buildindex)
end
