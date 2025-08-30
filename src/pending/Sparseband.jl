module Band

using ..MyData
#In v1.6 at sever, LinearAlgebra and BlockArrays are incompatible.
using LinearAlgebra
using BlockArrays

#Please use Energy function with Arpack very carefully, and use SparseArrays module carefully!!!
#Only eigs can handle sparse Matrices, and eigen cannot.
using Arpack
#similar work with sparseArray will return an "Array".
using SparseArrays


"""
`vector` = true or false(default), if use `vector` = true, recommend that run with MyWrite.WAVE.
"""

export Energy

function Energy(k::AbstractVector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real}; vector=true)

    H = Hamilton(k, hr, orbital, lattvec)

    if vector
        Ht = eigen(H)
        return MyEigen(Ht.values, Ht.vectors)
    else
        E = eigvals(H)
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end
function Hamilton(k::AbstractVector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real})
    norb = hr.norbital
    orblocat = orbital.location

    H = zeros(ComplexF64, norb, norb)

    for j = 1:norb
        for i = 1:j #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                dorb = orblocat[j, :] - orblocat[i, :]
                Ht = 0
                for ii in index
                    dR = dorb + hr.path[ii, 1] * lattvec[1, :] + hr.path[ii, 2] * lattvec[2, :] + hr.path[ii, 3] * lattvec[3, :]
                    Ht += hr.value[ii] * cis(k ⋅ dR)
                end
                H[i, j] = Ht
            end
        end
    end

    return Hermitian(H, :U)
end
################################################################################################################################################
function Energy(k::AbstractVector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real}, Nbands::Integer, fermienergy::Real; vector=true)

    H = Hamiltonsparse(k, hr, orbital, lattvec)

    Ht = eigs(H; nev=Nbands, sigma=fermienergy, ritzvec=vector)
    E = real(Ht[1])
    if maximum(abs, imag(Ht[1]) ./ E) > 1e-6
        println("Energy is not real!!!$(maximum(abs, imag(Ht[1]) ./ E))")
    end

    if vector
        return MyEigen(E, Ht[2])
    else
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end
function Hamilton(k::AbstractVector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real})
    norb = hr.norbital
    orblocat = orbital.location

    H = spzeros(ComplexF64, norb, norb)

    for j = 1:norb
        for i = 1:j #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                dorb = orblocat[j, :] - orblocat[i, :]
                Ht = 0
                for ii in index
                    dR = dorb + hr.path[ii, 1] * lattvec[1, :] + hr.path[ii, 2] * lattvec[2, :] + hr.path[ii, 3] * lattvec[3, :]
                    Ht += hr.value[ii] * cis(k ⋅ dR)
                end
                H[i, j] = Ht
            end
        end
    end

    return Hermitian(H, :U)
end
#############################################################################################################
function Energy(hr::HrData, orbital::OrbitalData; vector=true)

    H = Hamilton(hr, orbital)

    if vector
        Ht = eigen(H)
        return MyEigen(Ht.values, Ht.vectors)
    else
        E = eigvals(H)
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end

function Energy(hr::HrData, orbital::OrbitalData, Δ::Number; vector=true)

    H = HamiltonSC(hr, orbital, Δ)

    if vector
        Ht = eigen(H)
        return MyEigen(Ht.values, Ht.vectors)
    else
        E = eigvals(H)
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end

function Energy(hr::HrData, orbital::OrbitalData, Δ::Number, Nbands::Int, fermienergy::Real; vector=true)

    H = HamiltonSC(hr, orbital, Δ)

    Ht = eigs(H; nev=Nbands, sigma=fermienergy, ritzvec=vector)
    E = real(Ht[1])
    if maximum(abs, imag(Ht[1]) ./ E) > 1e-6
        println("Energy is not real!!!$(maximum(abs, imag(Ht[1]) ./ E))")
    end

    if vector
        return MyEigen(E, Ht[2])
    else
        return MyEigen(E, Matrix{ComplexF64}(undef, 0, 0))
    end
end

function HamiltonSC(hr::HrData, orbital::OrbitalData, Δ::Number)

    if Δ == 0
        H = Hamilton(hr, orbital)
    else

        He = Hamilton(hr, orbital) / 2
        Hh = -transpose(He)

        #Require the order of basis.
        norb = hr.norbital

        if issparse(He)

            iσy = sparse([0 Δ; -Δ 0]) / 2
            IΔ = blockdiagm(0 => repeat([iσy], norb ÷ 2))
            IΔ = sparse(IΔ)

            H = BlockArray(undef_blocks, SparseMatrixCSC{ComplexF64,Int}, [norb, norb], [norb, norb])
            H[Block(1, 1)] = He
            H[Block(2, 2)] = Hh
            H[Block(1, 2)] = IΔ
            H[Block(2, 1)] = IΔ'

            H = sparse(H)

        else

            iσy = [0 Δ; -Δ 0] / 2
            IΔ = blockdiagm(0 => repeat([iσy], norb ÷ 2))
            IΔ = Array(IΔ)

            H = BlockArray(undef_blocks, Matrix{ComplexF64}, [norb, norb], [norb, norb])
            H[Block(1, 1)] = He
            H[Block(2, 2)] = Hh
            H[Block(1, 2)] = IΔ
            H[Block(2, 1)] = IΔ'

            H = Array(H)

        end

    end

    return Hermitian(H, :U)
end

function Hamilton(hr::HrData, orbital::OrbitalData)
    norb = hr.norbital

    if issparse(hr.Nindex)
        H = spzeros(ComplexF64, norb, norb)
    else
        H = zeros(ComplexF64, norb, norb)
    end

    for j = 1:norb
        for i = 1:j #Accessed by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                Ht = 0
                for ii in index
                    Ht += hr.value[ii]
                end
                H[i, j] = Ht
            end
        end
    end

    return Hermitian(H, :U)
end

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

####################################################################################################
function Hamilton_all(k::Vector{Float64}, hr, orbital, lattvec::Matrix{Float64})
    norb = hr.norbital
    orblocat = orbital.location

    H = similar(hr.Nindex, ComplexF64)

    if isempty(hr.imagvalue)
        for j = 1:norb, i = 1:norb #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                dorb = orblocat[j, :] - orblocat[i, :]
                Ht = 0
                for ii in index
                    dR = dorb + hr.path[ii, 1] * lattvec[1, :] + hr.path[ii, 2] * lattvec[2, :] + hr.path[ii, 3] * lattvec[3, :]
                    Ht += hr.value[ii] * exp(im * (k ⋅ dR))
                end
                H[i, j] = Ht
            end
        end
    else
        for j = 1:norb, i = 1:norb #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                dorb = orblocat[j, :] - orblocat[i, :]
                Ht = 0
                for ii in index
                    dR = dorb + hr.path[ii, 1] * lattvec[1, :] + hr.path[ii, 2] * lattvec[2, :] + hr.path[ii, 3] * lattvec[3, :]
                    Ht += complex(hr.value[ii], hr.imagvalue[ii]) * exp(im * (k ⋅ dR))
                end
                H[i, j] = Ht
            end
        end
    end

    return (Hermitian(H, :U) + Hermitian(H, :L)) / 2
end

function Hamilton_all(hr, orbital, μ::Real)
    norb = hr.norbital

    H = H = similar(hr.Nindex, ComplexF64)

    if isempty(hr.imagvalue)
        for j = 1:norb, i = 1:norb #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                Ht = 0
                for ii in index
                    Ht += hr.value[ii]
                end
                H[i, j] = Ht
            end
        end
    else
        for j = 1:norb, i = 1:norb #Access by column.
            if hr.Nindex[i, j] > 0
                index = hr.index[i, j]
                Ht = 0
                for ii in index
                    Ht += complex(hr.value[ii], hr.imagvalue[ii])
                end
                H[i, j] = Ht
            end
        end
    end

    H = H - μ * I

    return (Hermitian(H, :U) + Hermitian(H, :L)) / 2
end

end
import .Band
