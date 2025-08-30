
function PlotBAND(kline::Kline, band::AbstractVector{<:Eigen}; kwards...)
	Nk = length(band)
	Nband = length(band[1].values)
	bandmatrix = Matrix{Float64}(undef, Nband, Nk)
	for i in eachindex(band)
		bandmatrix[:, i] .= band[i].values
	end
	return PlotBAND(kline, bandmatrix; kwards...)
end

function PlotBAND(kline::Kline, band::AbstractMatrix{<:Real};
	title = "",
	linewidth = 1.2,
	linecolor = :black,
	ylabel = "E/eV",
	xlims = kline.line[[1, end]] + [0, 0.001],
	ylims = (minimum(band), maximum(band)),
	xticks = (kline.line[kline.index], kline.name),
	legend = false,
	size = (600, 700),
	fermienergy = nothing,
)
	p = plot(kline.line, transpose(band);
		title,
		linewidth,
		linecolor,
		ylabel,
		xlims,
		ylims,
		xticks,
		legend,
		size,
	)
	if !isnothing(fermienergy)
		plot!(kline.line[[1, end]], fermienergy * ones(2); linecolor = :black, linestyle = :dot)
	end
	return p
end
