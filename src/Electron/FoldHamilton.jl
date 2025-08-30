export FoldHamilton, FoldUnitMat, FoldUnitTran!
function FoldHamilton(kdata::KData, Kₕ::KData, hr::HR, lattvec::AbstractMatrix{<:Real})

    norb = hr.norb
    Nk = kdata.N
    NKₕ = Kₕ.N

    Hk = Array{ComplexF64}(undef, norb * NKₕ, norb * NKₕ, Nk)

    for i in 1:Nk
        Hk[:, :, i] = blockdiagm(0 => map(k -> Array(Hamilton(k, hr, lattvec)), [kdata.point[i, :] + Kₕ.point[ii, :] for ii in 1:NKₕ]))
    end

    return Hk
end
#H(αₙ) = U'H(Gα)U
function FoldUnitMat(kmesh::KData, Kₕ::KData, lattvec::Matrix{<:Real}, norb::Integer)

    error("To be continued!")
    NK = Kₕ.N

    Kmeshindex = Kₕ.meshindex

    U = Matrix{Int}(undef, NK, NK)
    for i in 1:NK, j in 1:NK
        U[i, j] = Kmeshindex[i, :] ⋅ Kmeshindex[j, :]
    end

    dK = 2π / NK
    U = cis.(dK .* U) ./ sqrt(NK)

    I = diagm(ones(Int, norb))

    return kron(U, I)
end
function FoldUnitTran!(Hk::Array{<:Number}, kmesh::KData, Kₕ::KData, lattvec::Matrix{<:Real})

    NK = Kₕ.N

    Kmeshindex = Kₕ.meshindex

    U = Matrix{Int}(undef, NK, NK)
    for i in 1:NK, j in 1:NK
        U[i, j] = Kmeshindex[i, :] ⋅ Kmeshindex[j, :]
    end

    dK = 2π / NK
    U = cis.(dK .* U) ./ sqrt(NK)

    norb = size(Hk, 1) ÷ NK
    I = diagm(ones(Int, norb))

    U = kron(U, I)

    R = Kmeshindex * lattvec
    dR = Array{eltype(R)}(undef, NK, NK, 3)
    for i in 1:NK, j in 1:i
        dR[i, j, :] = R[i, :] - R[j, :]
    end


    I = ones(Int, norb, norb)
    Threads.@threads for i in 1:kmesh.N
        Uk = Matrix{ComplexF64}(undef, NK, NK)
        kpoint = kmesh.point[i, :]
        for ii in 1:NK, jj in 1:ii
            Uk[ii, jj] = cis(-kpoint ⋅ dR[ii, jj, :])
            Uk[jj, ii] = conj(Uk[ii, jj])
        end
        UUk = kron(Uk, I)
        Hk[:, :, i] = Hermitian((U' * Hk[:, :, i] * U) .* UUk, :U)
    end

    return Hk
end