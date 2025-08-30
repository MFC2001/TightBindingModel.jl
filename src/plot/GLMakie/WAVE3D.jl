export WAVE3D
function WAVE3D(location::AbstractVector, V::AbstractVector; kwards...)

	location_m = Matrix{Float64}(undef, 3, length(location))
	for (i, x) in enumerate(location)
		location_m[:, i] .= x
	end

	return WAVE3D(location_m', V; kwards...)
end
function WAVE3D(location::AbstractMatrix, V::AbstractVector;
	azimuth = -pi / 2,
	elevation = pi / 2,
)
	fig = GLMakie.Figure()
	ax = GLMakie.Axis3(
		fig[1, 1],
		aspect = :data,
		azimuth = azimuth,
		elevation = elevation,
		# xlabel = "",
		# ylabel = "",
		# zlabel = "",
		xticklabelsvisible = false,
		yticklabelsvisible = false,
		zticklabelsvisible = false,
		xticksvisible = false,
		yticksvisible = false,
		zticksvisible = false,
	)
	GLMakie.scatter!(ax,
		location[:, 1],
		location[:, 2],
		location[:, 3],
		color = :black,
		markersize = 5,
	)
	GLMakie.scatter!(ax,
		location[:, 1],
		location[:, 2],
		location[:, 3],
		color = :red,
		markersize = V,
	)
	return fig
end
