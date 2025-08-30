
#Pay attention to the kpoint where the band is degenerate.
function _check_symop_wave(irredkgrid::IrredKgrid, TB::AbstractTightBindModel; atol = 0.05)

	band_red = BAND(irredkgrid.redkdirect, TB; vector = true)
	band_irred = BAND(irredkgrid.kdirect, TB; vector = true)

	I = map(eachindex(irredkgrid.redkdirect)) do k
		ψ₁ = irredkgrid.lattsymop[k].U * band_irred[irredkgrid.irmap[k]].vectors
		ψ₂ = band_red[k].vectors
		all(_is_same_wave(ψ₁[:, i], ψ₂[:, i]; atol) for i in axes(ψ₁, 2))
	end

	return I
end
