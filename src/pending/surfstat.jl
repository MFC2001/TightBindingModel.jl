module SurfStatmodule
using LinearAlgebra
using ..MyData

export SurfStat
"""
    SurfStat(kpoints::AbstractMatrix{<:Real}, ω::AbstractVector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real}, direction::Integer;
    NP=1, outfolder::AbstractString="./", mode="surfstat", orbitalindex::AbstractVector{<:Pair}=["all" => 1:orbital.num])

`direction` = 1,2,3;
`mode` = "surfstat" or "greenfunc";
typeof(ω) <: AbstractVector{<:Real} or Real
"""
function SurfStat(kpoints::AbstractMatrix{<:Real}, ω::AbstractVector{<:Real}, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real}, direction::Integer;
    NP=1, η=abs(ω[1] - ω[end]) / length(ω), outfolder::AbstractString="./", mode="surfstat", orbitalindex::AbstractVector{<:Pair}=["all" => 1:orbital.num])

    if direction ∉ [1, 2, 3]
        error("Wrong direction from SurfStat.")
    end

    ωnum = length(ω)
    if ωnum == 1
        error("If only one ω value, please input ω instead of [ω].")
    end

    Nk = size(kpoints, 1)

    mkpath(outfolder)
    logfile = joinpath(outfolder, "surfstat.log")
    open(logfile, "w") do file
        println(file, "The number of kpoints is $(Nk).")
    end

    if mode == "surfstat"
        norb = length(orbitalindex)
        dos_l = Vector{Pair}(undef, norb)
        dos_r = Vector{Pair}(undef, norb)
        dos_bulk = Vector{Pair}(undef, norb)
        for i in 1:norb
            dos_l[i] = orbitalindex[i].first => Matrix{Float64}(undef, ωnum, Nk)
            dos_r[i] = orbitalindex[i].first => Matrix{Float64}(undef, ωnum, Nk)
            dos_bulk[i] = orbitalindex[i].first => Matrix{Float64}(undef, ωnum, Nk)
        end
    elseif mode == "greenfunc"
        norb = hr.norbital
        GLindex = 1:norb
        GRindex = (1:norb) .+ (NP - 1) * norb
    else
        error("Wrong mode from SurfStat.")
    end

    (hr, orbital, lattvec, orbitalindex_l, orbitalindex_r) = PrincipleLayer(hr, orbital, lattvec, NP, direction, orbitalindex)

    H00index = hr.path[:, direction] .== 0
    hr00 = HrData(hr.norbital, hr.path[H00index, :], hr.value[H00index])
    #-1 is corresponding to the function PrincipleLayer.
    H01index = hr.path[:, direction] .== -1
    path01 = hr.path[H01index, :]
    path01[:, direction] .= 0
    hr01 = HrData(hr.norbital, path01, hr.value[H01index])


    ω = ω .+ η * im

    if mode == "surfstat"
        for nk in axes(kpoints, 1)
            dt = @elapsed (GL_diag, GR_diag, GB_diag) = surfstat_kpoint(kpoints[nk, :], ω, hr00, hr01, orbital, lattvec)

            for i in 1:norb
                for j in 1:ωnum
                    dos_l[i].second[j, nk] = log(sum(GL_diag[orbitalindex_l[i].second, j]))
                    dos_r[i].second[j, nk] = log(sum(GR_diag[orbitalindex_r[i].second, j]))
                    dos_bulk[i].second[j, nk] = log(sum(GB_diag[orbitalindex_l[i].second, j]))
                end
            end

            open(logfile, "a") do file
                println(file, "kpoint $(nk) done, cost $(dt)s for irritation.")
            end
        end

        return dos_l, dos_r, dos_bulk

    elseif mode == "greenfunc"

        GL = Array{Float64}(undef, ωnum, Nk, norb)
        GR = Array{Float64}(undef, ωnum, Nk, norb)
        GB = Array{Float64}(undef, ωnum, Nk, norb)

        for nk in axes(kpoints, 1)
            dt = @elapsed (GL_diag, GR_diag, GB_diag) = surfstat_kpoint(kpoints[nk, :], ω, hr00, hr01, orbital, lattvec)

            GL[:, nk, :] = transpose(GL_diag[GLindex, :])
            GR[:, nk, :] = transpose(GR_diag[GRindex, :])
            GB[:, nk, :] = transpose(GB_diag[GLindex, :])

            open(logfile, "a") do file
                println(file, "kpoint $(nk) done, cost $(dt)s for irritation.")
            end
        end

        return GL, GR, GB
    end

end
function SurfStat(kpoints::AbstractMatrix{<:Real}, ω::Real, hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real}, direction::Integer;
    NP=1, η=0.002, outfolder::AbstractString="./", mode="surfstat", orbitalindex::AbstractVector{<:Pair}=["all" => 1:orbital.num])

    if direction ∉ [1, 2, 3]
        error("Wrong direction from SurfStat.")
    end

    Nk = size(kpoints, 1)

    mkpath(outfolder)
    logfile = joinpath(outfolder, "surfstat.log")
    open(logfile, "w") do
        println("The number of kpoints is $(Nk).")
    end

    if mode == "surfstat"
        norb = length(orbitalindex)
        dos_l = Vector{Pair}(undef, norb)
        dos_r = Vector{Pair}(undef, norb)
        dos_bulk = Vector{Pair}(undef, norb)
        for i in 1:norb
            dos_l[i] = orbitalindex[i].first => Matrix{Float64}(undef, 1, Nk)
            dos_r[i] = orbitalindex[i].first => Matrix{Float64}(undef, 1, Nk)
            dos_bulk[i] = orbitalindex[i].first => Matrix{Float64}(undef, 1, Nk)
        end
    elseif mode == "greenfunc"
        norb = hr.norbital
        GLindex = 1:norb
        GRindex = (1:norb) .+ (NP - 1) * norb
    else
        error("Wrong mode from SurfStat.")
    end

    (hr, orbital, lattvec, orbitalindex_l, orbitalindex_r) = PrincipleLayer(hr, orbital, lattvec, NP, direction, orbitalindex)

    H00index = hr.path[:, direction] .== 0
    hr00 = HrData(hr.norbital, hr.path[H00index, :], hr.value[H00index])
    #-1 is corresponding to the function PrincipleLayer.
    H01index = hr.path[:, direction] .== -1
    path01 = hr.path[H01index, :]
    path01[:, direction] .= 0
    hr01 = HrData(hr.norbital, path01, hr.value[H01index])


    ω = ω + η * im
    (GL_diag, GR_diag, GB_diag) = surfstat_ω(kpoints, ω, hr00, hr01, orbital, lattvec, logfile)

    if mode == "surfstat"
        for i in 1:norb
            for nk in axes(kpoints, 1)
                dos_l[i].second[1, nk] = log(sum(GL_diag[orbitalindex_l[i].second, nk]))
                dos_r[i].second[1, nk] = log(sum(GR_diag[orbitalindex_r[i].second, nk]))
                dos_bulk[i].second[1, nk] = log(sum(GB_diag[orbitalindex_l[i].second, nk]))
            end
        end

        return dos_l, dos_r, dos_bulk

    elseif mode == "greenfunc"
        GL = reshape(transpose(GL_diag[GLindex, :]), 1, Nk, norb)
        GR = reshape(transpose(GR_diag[GRindex, :]), 1, Nk, norb)
        GB = reshape(transpose(GB_diag[GLindex, :]), 1, Nk, norb)

        return GL, GR, GB
    end
end
###################################################################################################
function surfstat_kpoint(kpoint, ω, hr00, hr01, orbital, lattvec)
    H00 = layerH(kpoint, hr00, orbital, lattvec)
    H01 = layerH(kpoint, hr01, orbital, lattvec)

    ωnum = length(ω)
    norb = hr00.norbital
    GL_diag = Matrix{Float64}(undef, norb, ωnum)
    GR_diag = Matrix{Float64}(undef, norb, ωnum)
    GB_diag = Matrix{Float64}(undef, norb, ωnum)

    Threads.@threads for i in eachindex(ω)
        (GL, GR, GB) = surfgreen(ω[i], H00, H01)
        GL_diag[:, i] = -imag(diag(GL))
        GR_diag[:, i] = -imag(diag(GR))
        GB_diag[:, i] = -imag(diag(GB))
    end

    return GL_diag, GR_diag, GB_diag
end
function surfstat_ω(kpoints, ω, hr00, hr01, orbital, lattvec, logfile)

    Nk = size(kpoints, 1)
    norb = hr00.norbital
    GL_diag = Matrix{Float64}(undef, norb, Nk)
    GR_diag = Matrix{Float64}(undef, norb, Nk)
    GB_diag = Matrix{Float64}(undef, norb, Nk)

    Threads.@threads for i in axes(kpoints, 1)
        H00 = layerH(kpoints[i, :], hr00, orbital, lattvec)
        H01 = layerH(kpoints[i, :], hr01, orbital, lattvec)

        dt = @elapsed (GL, GR, GB) = surfgreen(ω, H00, H01)
        GL_diag[:, i] = -imag(diag(GL))
        GR_diag[:, i] = -imag(diag(GR))
        GB_diag[:, i] = -imag(diag(GB))

        open(logfile, "a") do file
            println(file, "kpoint $(i) done, cost $(dt)s for irritation.")
        end
    end

    return GL_diag, GR_diag, GB_diag
end
#####################################################################################################
function PrincipleLayer(hr::HrData, orbital::OrbitalData, lattvec::AbstractMatrix{<:Real}, NP::Integer, direction::Integer, orbitalindex)
    if NP == 1
        return hr, orbital, lattvec, orbitalindex, orbitalindex
    else
        norb = hr.norbital
        path = Matrix{Int}(undef, 0, 5)
        for i in 1:NP
            addpath = deepcopy(hr.path)
            addpath[:, 4] .+= norb * (i - 1)
            plrange = (1-NP:0) .+ (i - 1)
            layindex = map(x -> mod(x, plrange), addpath[:, direction])
            for j in 1:NP-1
                I = layindex .== plrange[j]
                addpath[I, 5] .+= norb * (NP - j)
            end
            addpath[:, direction] = (addpath[:, direction] - layindex) .÷ NP

            path = [path; addpath]
        end

        orbitalindex_bottom = similar(orbitalindex)
        for i in eachindex(orbitalindex)
            orbitalindex_bottom[i] = orbitalindex[i].first => (orbitalindex[i].second .+ norb * (NP - 1))
        end

        value = repeat(hr.value, NP)
        hr = HrData(norb * NP, path, value)

        orblocat = deepcopy(orbital.location)
        for i in 1:NP-1
            orblocat = [orblocat; orbital.location .- transpose(lattvec[direction, :] * i)]
        end
        orbital = OrbitalData(orblocat)

        lattvec = deepcopy(lattvec)
        lattvec[direction, :] .*= NP

        return hr, orbital, lattvec, orbitalindex, orbitalindex_bottom
    end
end
#####################################################################################################
function layerH(k, hr, orbital, lattvec)
    norb = hr.norbital
    orblocat = orbital.location

    H = zeros(ComplexF64, norb, norb)

    for j = 1:norb, i = 1:norb
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

    return H
end
##############################################################################################################
function surfgreen(ω::Complex, H00, H01; niter=101, err=1e-10)
    ϵl = H00
    ϵr = H00
    ϵ = H00
    α = H01
    β = H01'

    for i in 1:niter
        T = inv(ω * I - ϵ)
        Tα = α * T
        Tβ = β * T

        T = Tα * β
        ϵl = ϵl + T
        ϵ = ϵ + T

        T = Tβ * α
        ϵr = ϵr + T
        ϵ = ϵ + T

        α = Tα * α
        β = Tβ * β

        if norm(α) < err
            break
        elseif i == niter
            error("May can't converge.")
        end
    end

    GL = inv(ω * I - ϵl)
    GR = inv(ω * I - ϵr)
    GB = inv(ω * I - ϵ)

    return GL, GR, GB
end
end
using .SurfStatmodule