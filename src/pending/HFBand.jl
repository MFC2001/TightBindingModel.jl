module HFBand

using ..MyData

import SpecialFunctions: besselj
#In v1.6 at sever, LinearAlgebra and BlockArrays are incompatible.
using LinearAlgebra
using BlockArrays

#Please use Energy function with Arpack very carefully, and use SparseArrays module carefully!!!
#Only eigs can handle sparse Matrices, and eigen cannot.
using Arpack
#similar work with sparseArray will return an "Array".
using SparseArrays

#vector=true or false.
export Energy

# k -> k + qA/ħ
###############################################################################################################
struct Aparameter
    ampl::Real
    ω::Real # 10^15
    x::Real
    y::Real
    z::Real
    θx::Real
    θy::Real
end
#= In:
incident direction angle (θ, ϕ) in primary coordinate, acctually it's spherical coordinate of the opposite incident direction;
amplitude (Ax, Ay) in right-handed coordinate of light, the positive z-axis is opposite to the incident direction, and the x-axis is in the x-y plane of primary coordinate;
phase difference of incident light θ₀;
A = (Ax sin(ωt), Ay sin(ωt+θ₀)).
=#
#Out: A = (Ax sin(ωt+θx), Ay sin(ωt+θy), Az sin(ωt)) in primary coordinate.
function lightpara(θ::Real, ϕ::Real, Ax::Real, Ay::Real, ω::Real, θ₀::Real)::Aparameter
    Az = Ay * sin(θ)
    θz = θ₀

    #Ax
    Axx = -Ax * sin(ϕ) - Ay * cos(θ) * cos(ϕ) * cos(θ₀)
    Axy = -Ay * cos(θ) * cos(ϕ) * sin(θ₀)
    #Ay
    Ayx = Ax * cos(ϕ) - Ay * cos(θ) * sin(ϕ) * cos(θ₀)
    Ayy = -Ay * cos(θ) * sin(ϕ) * sin(θ₀)


    Ax = sqrt((Axx)^2 + (Axy)^2)
    θx = atan(Axy, Axx) - θz

    Ay = sqrt((Ayx)^2 + (Ayy)^2)
    θy = atan(Ayy, Ayx) - θz

    return Aparameter(max(Ax, Ay, Az), ω, Ax, Ay, Az, θx, θy)
end

##########################################################################################################################################
function Energy(k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real}, qmax::Integer, A::Aparameter; vector=true)

    if A.ampl < 1e-6
        H = Hamilton(k, hr, orbital, lattvec)
    else
        H = HF(k, hr, orbital, lattvec, qmax, A)
    end

    if vector
        Ht = eigen(H)

        norb = hr.norbital
        V = Ht.vectors[qmax*norb+1:(qmax+1)*norb, :]
        return MyEigen(Ht.values, V)
    else
        E = eigvals(H)
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end

function Energy(k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real}, qmax::Integer, A::Aparameter, Nbands::Integer, fermienergy::Real; vector=true)

    if A.ampl < 1e-6
        H = Hamilton(k, hr, orbital, lattvec)
    else
        H = HF(k, hr, orbital, lattvec, qmax, A)
    end

    Ht = eigs(H; nev=Nbands, sigma=fermienergy, ritzvec=vector)

    E = real(Ht[1])
    if maximum(abs, imag(Ht[1]) ./ E) > 1e-6
        error("Energy is not real!!!$(maximum(abs, imag(Ht[1]) ./ E))")
    end


    if vector
        norb = hr.norbital
        V = Ht[2][qmax*norb+1:(qmax+1)*norb, :]
        return MyEigen(E, V)
    else
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end

#########################################################################################################################
function Hamilton(k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real})
    norb = hr.norbital
    orblocat = orbital.location

    if issparse(hr.Nindex)
        H = spzeros(ComplexF64, norb, norb)
    else
        H = zeros(ComplexF64, norb, norb)
    end

    for j = 1:norb
        for i = 1:j #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                dorb = orblocat[j, :] .- orblocat[i, :]
                Ht = 0
                for ii in index
                    dR = dorb + hr.path[ii, 1] * lattvec[1, :] + hr.path[ii, 2] * lattvec[2, :]
                    Ht += hr.value[ii] * exp(im * (k ⋅ dR))
                end
                H[i, j] = Ht
            end
        end
    end

    return Hermitian(H, :U)
end

#################################################################################################################################
function HF(k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real}, qmax::Integer, A::Aparameter)
    if issparse(hr.Nindex)
        H = HFsparse(k, hr, orbital, lattvec, qmax, A)
    else
        H = HFarray(k, hr, orbital, lattvec, qmax, A)
    end
    return H
end

function HFarray(k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real}, qmax::Integer, A::Aparameter)
    norb = hr.norbital
    Nq = 2 * qmax + 1

    Hblock = zeros(ComplexF64, norb, norb)
    H = BlockArray(undef_blocks, Matrix{ComplexF64}, repeat([norb], Nq), repeat([norb], Nq))

    for n in -Nq+1:-qmax-1
        blockdiagm!(H, n => repeat([Hblock], Nq - abs(n)), -n => repeat([Hblock], Nq - abs(n)))
    end


    for nj in -qmax:-1
        Hnj = Hn(nj, k, hr, orbital, lattvec, A)
        # H_nj = Hn(-nj, k, hr, orblocat, lattvec, Apara)
        # Hnj = (Hnj + H_nj') / 2

        #= [   H₀    H_{nj}
            H_{-nj}    H₀  ] =#
        blockdiagm!(H, -nj => repeat([Hnj], Nq - abs(nj)), nj => repeat([Hnj'], Nq - abs(nj)))
    end

    #[-qmax, -qmax+1, ... 0, ... qmax]ᵀ.
    H₀ = Hn(0, k, hr, orbital, lattvec, A)
    for i in -qmax:qmax
        n = i + qmax + 1
        H[Block(n, n)] = H₀ + i * A.ω * I
    end

    return Hermitian(Array(H), :U)
end

function HFsparse(k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real}, qmax::Integer, A::Aparameter)
    norb = hr.norbital
    Nq = 2 * qmax + 1

    Hblock = spzeros(ComplexF64, norb, norb)
    H = BlockArray(undef_blocks, SparseMatrixCSC{ComplexF64,Int}, repeat([norb], Nq), repeat([norb], Nq))

    for n in -Nq+1:-qmax-1
        blockdiagm!(H, n => repeat([Hblock], Nq - abs(n)), -n => repeat([Hblock], Nq - abs(n)))
    end

    for nj in -qmax:-1
        Hnj = Hn(nj, k, hr, orbital, lattvec, A)
        # H_nj = Hn(-nj, k, hr, orblocat, lattvec, Apara)
        # Hnj = (Hnj + H_nj') / 2

        #= [   H₀    H_{nj}
            H_{-nj}    H₀  ] =#
        blockdiagm!(H, -nj => repeat([Hnj], Nq - abs(nj)), nj => repeat([Hnj'], Nq - abs(nj)))
    end

    #[-qmax, -qmax+1, ... 0, ... qmax]ᵀ.
    H₀ = Hn(0, k, hr, orbital, lattvec, A)
    for i in -qmax:qmax
        n = i + qmax + 1
        H[Block(n, n)] = H₀ + i * A.ω * I
    end

    return Hermitian(sparse(H), :U)
end

# k -> k + qA/ħ
#qA/ħ = (Ax sin(ωt+θx), Ay sin(ωt+θy), Az sin(ωt)).
function Hn(nj::Integer, k::Vector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::Matrix{<:Real}, A::Aparameter)
    orblocat = orbital.location

    αr(rx, ry, rz) = A.x * rx * cos(A.θx) + A.y * ry * cos(A.θy) + A.z * rz
    βr(rx, ry, rz) = A.x * rx * sin(A.θx) + A.y * ry * sin(A.θy)
    function h(r)
        α = αr(r...)
        β = βr(r...)
        return besselj(-nj, sqrt(α^2 + β^2)) * cis(k ⋅ r - nj * atan(β, α))
    end

    norb = hr.norbital
    if issparse(hr.Nindex)
        H = spzeros(ComplexF64, norb, norb)
    else
        H = zeros(ComplexF64, norb, norb)
    end

    for j = 1:norb, i = 1:norb #Access by column.
        if hr.Nindex[i, j] > 0
            index = hr.index[i, j]
            dorb = orblocat[j, :] .- orblocat[i, :]
            Ht = 0
            for a in index
                dR = dorb + hr.path[a, 1] * lattvec[1, :] + hr.path[a, 2] * lattvec[2, :]
                Ht += hr.value[a] * h(dR)
            end
            H[i, j] = Ht
        end
    end

    return H
end

##########################################################################################################
function blockdiagm(args::Pair{<:Integer,<:AbstractVector}...)
    N = length(args[1].second) + abs(args[1].first)

    blockdim = ndims(args[1].second[1])
    blocksize = size(args[1].second[1])
    if blockdim == 1
        a = blocksize[1]
        b = 1
        blocktype = typeof(similar(args[1].second[1], 1, 1))

        blockmatrix = BlockArray(undef_blocks, blocktype, repeat([a], N), repeat([b], N))

        index = Vector{Integer}(undef, 0)
        for arg in args
            n = arg.first
            blockindex = diagm_Blockindex(N, n)
            for i in eachindex(blockindex)
                blockmatrix[blockindex[i]] = reshape(arg.second[i], :, 1)
            end
            push!(index, n)
        end

    elseif blockdim == 2
        a = blocksize[1]
        b = blocksize[2]
        blocktype = typeof(args[1].second[1])

        blockmatrix = BlockArray(undef_blocks, blocktype, repeat([a], N), repeat([b], N))

        index = Vector{Integer}(undef, 0)
        for arg in args
            n = arg.first
            blockindex = diagm_Blockindex(N, n)
            for i in eachindex(blockindex)
                blockmatrix[blockindex[i]] = arg.second[i]
            end
            push!(index, n)
        end
    else
        error("Wrong blockdim, only work for Vector or Matrix block!")
    end

    if issparse(args[1].second[1])
        T = spzeros(eltype(args[1].second[1]), a, b)
    else
        T = zeros(eltype(args[1].second[1]), a, b)
    end

    for n in setdiff(-N+1:N-1, index)
        for blockindex in diagm_Blockindex(N, n)
            blockmatrix[blockindex] = T
        end
    end

    return blockmatrix
end
function blockdiagm!(blockmatrix::AbstractBlockMatrix, args::Pair{<:Integer,<:AbstractVector}...)
    N = length(args[1].second) + abs(args[1].first)

    blockdim = ndims(args[1].second[1])
    blocksize = size(args[1].second[1])
    if blockdim == 1
        a = blocksize[1]
        b = 1
        blocktype = typeof(similar(args[1].second[1], 1, 1))

        index = Vector{Integer}(undef, 0)
        for arg in args
            n = arg.first
            blockindex = diagm_Blockindex(N, n)
            for i in eachindex(blockindex)
                blockmatrix[blockindex[i]] = reshape(arg.second[i], :, 1)
            end
            push!(index, n)
        end

    elseif blockdim == 2
        a = blocksize[1]
        b = blocksize[2]
        blocktype = typeof(args[1].second[1])

        index = Vector{Integer}(undef, 0)
        for arg in args
            n = arg.first
            blockindex = diagm_Blockindex(N, n)
            for i in eachindex(blockindex)
                blockmatrix[blockindex[i]] = arg.second[i]
            end
            push!(index, n)
        end
    else
        error("Wrong blockdim, only work for Vector or Matrix block!")
    end

    return nothing
end
function diagm_Blockindex(N::Integer, n::Integer)
    index = Vector{Block{2,Int64}}(undef, 0)

    if n >= 0
        for i in 1:N-abs(n)
            push!(index, Block(i, i + n))
        end
    else
        for i in 1:N-abs(n)
            push!(index, Block(i - n, i))
        end
    end

    return index
end

end
import .HFBand