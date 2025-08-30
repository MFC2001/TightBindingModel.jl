export V2DRK
"""
2D Rytova–Keldysh Potential is another version of a Keldysh type Coulomb potential used for 2D systems.
https://doi.org/10.1103/physrevb.98.125308
"""
function V2DRK(lattvec::AbstractMatrix{<:Real}; LC=nothing, w₀=0, r₀=1, ϵ=1, ϵz=1, ϵ₁=1, ϵ₂=1, q_tol=0.001)::Function
    #unit: Å, 2D model is in x-y plane.
    # w₀ and r₀ are free parameters;
    # ϵ are in-plane macroscopic effective dielectric constant, ϵz is the out of plane effective dielectric constant of the crystal;
    # ϵ₁ and ϵ₂ correspond to the dielectric constants of the substrate above and below.

    #1e-19
    qₑ = 1.602176634
    #1e-12
    ϵ₀ = 8.854187817

    η = √(ϵ / ϵz)
    κ = √(ϵ * ϵz)

    p₁ = (ϵ₁ - κ) / (ϵ₁ + κ)
    p₂ = (ϵ₂ - κ) / (ϵ₂ + κ)


    a₁ = lattvec[:, 1]
    a₂ = lattvec[:, 2]
    Auc = norm(a₁ × a₂)

    if isnothing(LC)
        LC = lattvec[3, 3]
    end


    #E/eV, so use qₑ istead of qₑ^2.
    T = qₑ * 1e3 / (2 * Auc * ϵ₀)

    F(q) = (1 - p₁ * p₂ * exp(-2 * q * η * LC)) * κ * exp(w₀ * q) / (1 - p₁ * exp(-η * q * LC)) / (1 - p₂ * exp(-η * q * LC)) + r₀ * q
    V = function (q)
        q = norm(q)
        return q < q_tol ? 0 : T / (q * F(q))
    end

    return V
end