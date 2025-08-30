module Exciton

using LinearAlgebra
using ..MyData
using ..BAND
import ..MyFBZ: kmeshmap

using NLsolve

function main(kmesh::KData, hr::HrData, orbital::OrbitalData, POSCAR::POSCARData; outfolder::AbstractString="./", MV::Integer, Δeps=1e-10, μ=0, T=0, lc=nothing, ϵ=1)

    f = FDdistribution(μ, T)

    Nk = kmesh.N

    original_band = Band(kmesh, hr, POSCAR.lattice_vec)

    V = Coulomb(; ϵ)
    Vq = Vmatrix(V, POSCAR, orbital, kmesh)

    addmap = kmeshmap(kmesh; mode=+)
    minusmap = kmeshmap(kmesh; mode=-)

    K = kernal(Vq, original_band, MV; addmap, minusmap)

    norb = hr.norbital

    band = Vector{Eigen}(undef, Nk)
    Threads.@threads for i in 1:Nk
        H = diagm(complex.(original_band[i].values))
        H = Hermitian(H, :U)
        band[i] = eigen!(H)
    end
    Δ₀ = orderpara(f, K, band; addmap, minusmap)

    # Threads.@threads for i in 1:Nk
    #     H = diagm(complex.(original_band[i].values)) + Δ₀[:, :, i]
    #     H = diagm(complex.(original_band[i].values))
    #     H[1:MV, MV+1:norb] = rand(Float64, MV, MV) .% 0.2
    #     # H[1:MV, MV+1:norb] .= 0.1
    #     H = Hermitian(H, :U)
    #     band[i] = eigen!(H)
    # end
    # Δ₀ = orderpara(f, K, band; addmap, minusmap)

    # F! = function (out, x)
    #     Δ = x2Δ(x, norb, Nk)
    #     Threads.@threads for i in 1:Nk
    #         H = Hermitian(diagm(complex.(original_band[i].values)) + Δ[:, :, i], :U)
    #         band[i] = eigen!(H)
    #     end

    #     Δ = orderpara(f, K, band; addmap, minusmap)

    #     out .= x - Δ2x(Δ)
    #     return nothing
    # end

    # @time result = nlsolve(F!, Δ2x(Δ₀), method=:newton, ftol=Δeps, show_trace=true)
    # Δ = x2Δ(result.zero, norb, Nk)
    Δ = Δ₀

    Hk = Array{ComplexF64}(undef, norb, norb, Nk)
    for i in 1:Nk
        H = Hermitian(diagm(complex.(original_band[i].values)) + Δ[:, :, i], :U)
        U = original_band[i].vectors
        H = U * H * U'
        Hk[:, :, i] = H
    end

    hr = inversefourier(Hk, kmesh, POSCAR.lattice_vec, 1e-6)

    return hr
end

###########################################################################################
function Δ2x(Δ)
    norb = size(Δ, 1)
    Nk = size(Δ, 3)
    x = Vector{real(eltype(Δ))}(undef, Nk * norb * norb)

    t = 1
    for k in 1:Nk
        for i in 2:norb
            for j in 1:i-1
                d = Δ[i, j, k]
                x[t] = real(d)
                x[t+1] = imag(d)
                t += 2
            end
        end
        for i in 1:norb
            x[t] = real(Δ[i, i, k])
            t += 1
        end
    end

    return x
end
function x2Δ(x, norb, Nk)
    Δ = Array{complex(eltype(x))}(undef, norb, norb, Nk)
    TΔ = Matrix{complex(eltype(x))}(undef, norb, norb)

    n = norb * norb
    for k in 1:Nk
        T = @view x[n*(k-1)+1:n*k]
        t = 1
        for i in 2:norb
            for j in 1:i-1
                TΔ[i, j] = complex(T[t], T[t+1])
                t += 2
            end
        end
        for i in 1:norb
            TΔ[i, i] = complex(T[t])
            t += 1
        end
        Δ[:, :, k] = Hermitian(TΔ, :L)
    end

    return Δ
end

#######################################################################
function orderpara(f, K, band; addmap, minusmap)
    nv = size(K, 1)
    nc = size(K, 2)
    norb = length(band[1].values)
    Nk = length(band)


    DM = DensityMatrix(f, band)

    Δ = Array{ComplexF64}(undef, norb, norb, Nk)
    Threads.@threads for k in 1:Nk
        TΔ = zeros(ComplexF64, norb, norb)
        # for v in 1:nv, c in 1:nc
        #     for i in 1:Nk
        #         TΔ[v, nv+c] += sum(K[:, :, i, v, c, k, minusmap[1, 1]] .* transpose(DM[nv+1:end, 1:nv, i]))
        #     end
        # end

        # for v1 in 1:nv, v2 in v1:nv
        #     for q in 1:Nk
        #         TΔ[v1, v2] -= sum(K[v2, :, k, v1, :, k, q] .* DM[nv+1:end, nv+1:end, addmap[k, q]])
        #     end
        # end
        # for c1 in 1:nc, c2 in c1:nc
        #     for q in 1:Nk
        #         kq = minusmap[k, q]
        #         TΔ[nv+c1, nv+c2] -= sum(K[:, c1, kq, :, c2, kq, q] .* transpose(DM[1:nv, 1:nv, kq]))
        #     end
        # end

        for v in 1:nv, c in 1:nc
            for i in 1:Nk
                TΔ[v, nv+c] += sum(K[:, :, i, v, c, k] .* transpose(DM[nv+1:end, 1:nv, i]))
            end
        end

        for v1 in 1:nv, v2 in v1:nv
            TΔ[v1, v2] -= sum(K[v2, :, k, v1, :, k] .* DM[nv+1:end, nv+1:end, k])
        end
        for c1 in 1:nc, c2 in c1:nc
            TΔ[nv+c1, nv+c2] -= sum(K[:, c1, k, :, c2, k] .* transpose(DM[1:nv, 1:nv, k]))
        end

        Δ[:, :, k] = Hermitian(TΔ, :U)
    end

    return Δ
end
##################################################################
function DensityMatrix(f, band)
    norb = length(band[1].values)
    Nk = length(band)

    DM = Array{ComplexF64}(undef, norb, norb, Nk)
    Threads.@threads for k in 1:Nk
        fvector = f.(band[k].values)
        wavek = band[k].vectors
        fwave = wavek .* transpose(fvector)
        for i in 2:norb, j in 1:i-1
            DM[i, j, k] = wavek[i, :] ⋅ fwave[j, :]
            DM[j, i, k] = conj(DM[i, j, k])
        end
        for i in 1:norb
            DM[i, i, k] = real(wavek[i, :] ⋅ fwave[i, :])
        end
    end

    return DM
end
##########################################################################
function kernal(V, original_band, MV; addmap, minusmap)
    nband = length(original_band[1].values)
    nv = MV
    nc = nband - nv
    Nk = length(original_band)

    # T1 = Array{ComplexF64}(undef, nband, nband, nv, nc, Nk, Nk, Nk)
    # Threads.@threads for q in 1:Nk
    #     for i in 1:nband, j in 1:nband, v in 1:nv, c in 1:nc, k1 in 1:Nk, k2 in 1:Nk
    #         T1[i, j, v, c, k1, k2, q] = original_band[k1].vectors[i, v] * original_band[addmap[k2, q]].vectors[j, c+nv]
    #     end
    # end

    # #K v'c'k' vck q
    # K = Array{ComplexF64}(undef, nv, nc, Nk, nv, nc, Nk, Nk)
    # Threads.@threads for q in 1:Nk
    #     for v1 in 1:nv, c1 in 1:nc, k1 in 1:Nk, v2 in 1:nv, c2 in 1:nc, k2 in 1:Nk
    #         t = 0
    #         for i in 1:nband, j in 1:nband
    #             t += conj(T1[i, j, v2, c1, k2, k1, q]) * (T1[j, i, v1, c2, k1, k2, q] * V[i, j, minusmap[minusmap[1, 1], q]]
    #                                                       -
    #                                                       T1[i, j, v1, c2, k1, k2, q] * V[i, j, minusmap[k2, k1]])
    #             # t -= conj(T1[i, j, v2, c1, k2, k1, q]) * T1[i, j, v1, c2, k1, k2, q] * V[i, j, minusmap[k2, k1]]
    #         end
    #         K[v1, c1, k1, v2, c2, k2, q] = t
    #     end
    # end

    T1 = Array{ComplexF64}(undef, nband, nband, nv, nc, Nk, Nk)
    Threads.@threads for k2 in 1:Nk
        for i in 1:nband, j in 1:nband, v in 1:nv, c in 1:nc, k1 in 1:Nk
            T1[i, j, v, c, k1, k2] = original_band[k1].vectors[i, v] * original_band[k2].vectors[j, c+nv]
        end
    end

    #K v'c'k' vck
    K = Array{ComplexF64}(undef, nv, nc, Nk, nv, nc, Nk)
    Threads.@threads for k2 in 1:Nk
        for v1 in 1:nv, c1 in 1:nc, k1 in 1:Nk, v2 in 1:nv, c2 in 1:nc
            t = 0
            for i in 1:nband, j in 1:nband
                t += conj(T1[i, j, v2, c1, k2, k1]) * (T1[j, i, v1, c2, k1, k2] * V[i, j, minusmap[1, 1]]
                                                       -
                                                       T1[i, j, v1, c2, k1, k2] * V[i, j, minusmap[k2, k1]])
                # t -= conj(T1[i, j, v2, c1, k2, k1]) * T1[i, j, v1, c2, k1, k2] * V[i, j, minusmap[k2, k1]]
            end
            K[v1, c1, k1, v2, c2, k2] = t
        end
    end

    return K
end

#############################################################################
function Vmatrix(V, POSCAR, orbital, kmesh)
    norb = orbital.num
    kpoint = kmesh.point
    Nk = kmesh.N

    ucpath = getucpath(Nk, POSCAR.periodicity)
    N = size(ucpath, 1)
    R = ucpath * POSCAR.lattice_vec

    orblocat = orbital.location
    B = Array{Float64}(undef, norb, norb, 3)
    for i in 1:norb, j in 1:i
        dorb = orblocat[i, :] - orblocat[j, :]
        B[i, j, :] = dorb
        B[j, i, :] = -dorb
    end

    VR = Array{Float64}(undef, norb, norb, N)
    Threads.@threads for i in 1:N
        R₀ = R[i, :]
        for ii in 1:norb, jj in 1:norb
            VR[ii, jj, i] = V(R₀ + B[ii, jj, :])
        end
    end

    Vq = zeros(ComplexF64, norb, norb, Nk)
    Threads.@threads for i in 1:Nk
        k = kpoint[i, :]
        for ii in 1:norb, jj in 1:norb
            T = 0
            for j in 1:N
                T += VR[ii, jj, j] * cis(k ⋅ R[j, :])
            end
            Vq[ii, jj, i] = T / Nk
        end
    end

    return Vq
end
function getucpath(Nk, periodicity)

    p = count(periodicity .== "p")

    if p == 3
        shell = Int(ceil(∛Nk / 2))

        T = [repeat(-shell:shell, inner=2 * shell + 1) repeat(-shell:shell, outer=2 * shell + 1)]
        ucpath = [repeat(T, inner=(2 * shell + 1, 1)) repeat(-shell:shell, outer=(2 * shell + 1)^2)]

    elseif p == 2
        circle = Int(ceil(√Nk / 2))
        ucpath = zeros(Int, (2 * circle + 1)^2, 3)

        if periodicity[1] == "np"
            T = 2:3
        elseif periodicity[2] == "np"
            T = [3, 1]
        elseif periodicity[3] == "np"
            T = 1:2
        else
            error("Wrong periodicity from ExpandSuperCell!.")
        end

        ucpath[:, T] = [repeat(-circle:circle, inner=2 * circle + 1) repeat(-circle:circle, outer=2 * circle + 1)]

    elseif p == 1
        t = Int(ceil(Nk / 2))
        ucpath = zeros(Int, 2 * circle + 1, 3)

        if periodicity[1] == "p"
            T = 1
        elseif periodicity[2] == "p"
            T = 2
        elseif periodicity[3] == "p"
            T = 3
        else
            error("Wrong periodicity from inversefourier.")
        end
        ucpath[:, T] = -t:t
    elseif p == 0
        ucpath = [0 0 0]
    else
        error("Wrong periodicity from inversefourier.")
    end

    return ucpath
end

function Coulomb(; ϵ=1)

    #1e-19
    qₑ = 1.602176634
    #1e-12
    ϵ₀ = 8.854187817

    T = qₑ * 1e3 / (4 * π * ϵ * ϵ₀)

    V = function (r)
        r = norm(r)
        return r < 3 ? T / 3 : T / norm(r)
    end

    return V
end

#2D Keldysh Potential(V2DK)
function V2DK(Nk::Integer, POSCAR::POSCARData; lc=nothing, ϵ₁=1, ϵ₂=1, ϵ₃=1, q_tol=0.001)::Function
    #lc/A

    #1e-19
    qₑ = 1.602176634
    #1e-12
    ϵ₀ = 8.854187817
    ϵd = (ϵ₁ + ϵ₃) / 2


    lattvec = POSCAR.lattice_vec
    periodicity = POSCAR.periodicity

    p = count(periodicity .== "p")
    if p == 2
        if periodicity[1] == "np"
            a₁ = lattvec[2, :]
            a₂ = lattvec[3, :]
            a₃ = lattvec[1, :]
        elseif periodicity[2] == "np"
            a₁ = lattvec[1, :]
            a₂ = lattvec[3, :]
            a₃ = lattvec[2, :]
        elseif periodicity[3] == "np"
            a₁ = lattvec[1, :]
            a₂ = lattvec[2, :]
            a₃ = lattvec[3, :]
        else
            error("Wrong periodicity from kline.")
        end
    else
        error("Wrong periodicity, it should have 2 \"p\".")
    end

    T = a₁ × a₂
    Auc = norm(T)
    a = (norm(a₁) + norm(a₂)) / 2
    if isnothing(lc)
        lc = (a₃ ⋅ T) / Auc
    end

    r₀ = lc * (ϵ₂ - 1) / (ϵ₁ + ϵ₃)

    Δ = 2 * π * r₀ / (a * Nk)

    α₁ = 1.76
    α₂ = 1
    α₃ = 0

    #E/eV, so use qₑ istead of qₑ^2.
    T = qₑ * 1e3 / (2 * Nk * Auc * ϵ₀ * ϵd)
    T1 = qₑ * a * 1e3 / (4π * Auc * ϵ₀ * ϵd)

    W(q) = begin
        q = norm(q)
        q < q_tol ? T1 * (α₁ + α₂ * Δ + α₃ * Δ^2) : T / (q * (1 + r₀ * q))
    end

    return W
end

#####################################################################################
function FDdistribution(μ::Real, T::Real)::Function
    #energy unit is eV, so use kB/e.
    kB = 8.6173332621e-5

    if T == 0
        f0(ϵ) = ϵ > μ ? 0 : 1
        return f0
    elseif T > 0
        f(ϵ) = 1 / (exp((ϵ - μ) / (kB * T)) + 1)
        return f
    else
        error("Tempreture should be not less than 0. ")
    end
end

end
import .Exciton