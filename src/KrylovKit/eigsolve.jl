
function _eigsolve_Hmat(H)

	(vals, vecs, info) = eigsolve(H, eigsolveconfig[:howmany], eigsolveconfig[:which];
		verbosity = eigsolveconfig[:verbosity],
		tol = eigsolveconfig[:tol],
		krylovdim = min(eigsolveconfig[:krylovdim], size(H, 1)),
		maxiter = eigsolveconfig[:maxiter],
		orth = eigsolveconfig[:orth],
		ishermitian = true,
	)

	if info.converged < eigsolveconfig[:howmany]
		@warn "KrylovKit.eigsolve don't converge!"
	elseif info.converged > eigsolveconfig[:howmany]
		vals = vals[1:eigsolveconfig[:howmany]]
		vecs = vecs[1:eigsolveconfig[:howmany]]
	end

	I = sortperm(vals)
	vals = vals[I]

	vecs_mat = Matrix{eltype(eltype(vecs))}(undef, size(H, 1), eigsolveconfig[:howmany])
	for i in eachindex(I)
		vecs_mat[:, i] = vecs[I[i]]
	end

	return Eigen(vals, vecs_mat)
end

function _eigen2vals(eigens::AbstractArray{<:Eigen})
	isempty(eigens) && return Matrix{Int}(undef, 0, 0)
	len_values = [length(e.values) for e in eigens]
	max_len = maximum(len_values)
	n = length(eigens)
	value = fill(convert(eltype(first(eigens).values), NaN), max_len, n)
	# Matrix{eltype(first(eigens).values)}(undef, max_len, n)
	for (i, e) in enumerate(eigens)
		value[1:len_values[i], i] .= e.values
	end
	return value
end
