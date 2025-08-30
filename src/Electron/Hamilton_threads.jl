export Hamilton, HamiltonSC
function Hamilton(k::AbstractVector{<:Real}, hr::HR, lattvec::AbstractMatrix{<:Real})
    kl = transpose(lattvec) * k

    ittr = Vector{Tuple{Int,Int}}(undef, Int(hr.norb * (hr.norb + 1) / 2))
    n = 0
    for j in 1:hr.norb, i in 1:j
        n += 1
        ittr[n] = (i, j)
    end

    H = Matrix{ComplexF64}(undef, hr.norb, hr.norb)

    Threads.@threads for (i, j) in ittr
        if hr.Nindex[i, j] > 0
            H[i, j] = sum(hr.index[i, j]) do index
                return hr.value[index] * cis(kl ⋅ hr.path[index])
            end
        else
            H[i, j] = 0
        end
    end

    return Hermitian(H, :U)
end
function Hamilton(k::AbstractVector{<:Real}, hr::HR, orbital::ORBITAL, lattvec::AbstractMatrix{<:Real})
    kl = transpose(lattvec) * k
    orblocat = orbital.location

    ittr = Vector{Tuple{Int,Int}}(undef, Int(hr.norb * (hr.norb + 1) / 2))
    n = 0
    for j in 1:hr.norb, i in 1:j
        n += 1
        ittr[n] = (i, j)
    end

    H = Matrix{ComplexF64}(undef, hr.norb, hr.norb)

    Threads.@threads for (i, j) in ittr
        if hr.Nindex[i, j] > 0
            dorb = orblocat[j] - orblocat[i]
            H[i, j] = cis(k ⋅ dorb) * sum(hr.index[i, j]) do index
                return hr.value[index] * cis(kl ⋅ hr.path[index])
            end
        else
            H[i, j] = 0
        end
    end

    return Hermitian(H, :U)
end

#############################################################################
function HamiltonSC(hr::HR, Δ::Number)

    if Δ == 0
        H = Hamilton(hr)
    else

        He = Hamilton(hr) / 2
        Hh = -transpose(He)

        #Require the order of basis.
        norb = hr.norb


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

    return Hermitian(H, :U)
end

function Hamilton(hr::HR)

    ittr = Vector{Tuple{Int,Int}}(undef, Int(hr.norb * (hr.norb + 1) / 2))
    n = 0
    for j in 1:hr.norb, i in 1:j
        n += 1
        ittr[n] = (i, j)
    end

    H = Matrix{ComplexF64}(undef, hr.norb, hr.norb)

    Threads.@threads for (i, j) in ittr
        if hr.Nindex[i, j] > 0
            H[i, j] = sum(hr.value[hr.index[i, j]])
        else
            H[i, j] = 0
        end
    end

    return Hermitian(H, :U)
end
