export V2DT2
"""
A variation for the truncated Coulomb potential for 2D sys-tems proposed by IsmailBeigi et al. https://doi.org/10.1103/physrevb.73.233103
"""
function V2DT2(lattvec::AbstractMatrix{<:Real}; LC=nothing, q_tol=0.001)::Function
    #unit: Å, 2D model is in x-y plane.

    Vuc = abs((lattvec[:, 1] × lattvec[:, 2]) ⋅ lattvec[:, 3])

    if isnothing(LC)
        ZC = lattvec[3, 3] / 2
    else
        ZC = LC / 2
    end

    #1e-19
    qₑ = 1.602176634
    #1e-12
    ϵ₀ = 8.854187817

    T = qₑ * 1e3 / (2 * Vuc * ϵ₀)
    q_tol2 = q_tol^2

    V = function (Q)
        q2 = sum(abs2, Q)
        return q2 < q_tol2 ? 0 : T / q2 * (1 - exp(-√(q2 - Q[3]^2) * ZC) * cos(Q[3] * ZC))
    end

    return V
end