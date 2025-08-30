export V2DOH
"""
2D Rytova–Keldysh Potential is another version of a Keldysh type Coulomb potential used for 2D systems.
https://doi.org/10.1103/physrevb.98.125308
"""
function V2DOH(lattvec::AbstractMatrix{<:Real}; w₀=0, ϵ=1, q_tol=0.001)::Function
    #unit: Å, 2D model is in x-y plane.
    # w₀ is free parameter;
    # ϵ is an effective dielectric constant parameter;

    #1e-19
    qₑ = 1.602176634
    #1e-12
    ϵ₀ = 8.854187817


    a₁ = lattvec[:, 1]
    a₂ = lattvec[:, 2]
    Auc = norm(a₁ × a₂)


    #E/eV, so use qₑ istead of qₑ^2.
    T = qₑ * 1e3 / (2 * Auc * ϵ₀ * ϵ)

    V = function (q)
        q = norm(q)
        return q < q_tol ? 0 : T / (q * exp(w₀ * q))
    end

    return V
end