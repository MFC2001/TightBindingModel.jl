export V2DK
"""
2D Keldysh Potential(V2DK). https://doi.org/10.1103/physrevb.97.205409
"""
function V2DK(lattvec::Lattice, Nk::Integer; LC=nothing, ϵ=1, ϵ₁=1, ϵ₂=1, q_tol=0.001)::Function
    #unit: Å, 2D model is in x-y plane.
    # ϵ are in-plane macroscopic effective dielectric constant;
    # ϵ₁ and ϵ₂ correspond to the dielectric constants of the substrate above and below.

    #1e-19
    qₑ = 1.602176634
    #1e-12
    ϵ₀ = 8.854187817
    ϵd = (ϵ₁ + ϵ₂) / 2

    #https://doi.org/10.1103/physrevb.97.205409
    α₁ = 1.76
    α₂ = 1
    α₃ = 0


    a₁ = lattvec[:, 1]
    a₂ = lattvec[:, 2]

    Auc = norm(a₁ × a₂)
    a = (norm(a₁) + norm(a₂)) / 2

    if isnothing(LC)
        LC = lattvec[3, 3]
    end

    r₀ = LC * (ϵ - 1) / (ϵ₁ + ϵ₂)
    Δ = 2 * π * r₀ / (a * √Nk)



    #E/eV, so use qₑ istead of qₑ^2.
    T = qₑ * 1e3 / (2 * Auc * ϵ₀ * ϵd)
    V₀ = qₑ * a * √Nk * 1e3 / (4π * Auc * ϵ₀ * ϵd) * (α₁ + α₂ * Δ + α₃ * Δ^2)

    V = function (q)
        q = norm(q)
        return q < q_tol ? V₀ : T / (q * (1 + r₀ * q))
    end

    return V
end