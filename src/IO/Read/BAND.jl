export ReadBAND, smoothband, wannierband
export WAVE

"""
    ReadBAND(file::AbstractString; readkname='N', smoothband='N')

Read BAND.dat.
"""
function ReadBAND(file::AbstractString; readkname='N', smoothband='N')

    file = open(file, "r")
    readline(file)

    line = readline(file)
    line = split(line, ':')[2]
    line = parse.(Int, split(line))
    Nk = line[1]
    Nband = line[2]

    if readkname[1] ∈ ['N', 'n']
        kindex = Int[]
        kname = String[]
    elseif readkname[1] ∈ ['Y', 'y']
        #kindex.
        line = readline(file)
        line = split(line, ':')[2]
        kindex = parse.(Int, split(line))
        #kpoint
        readline(file)
        #kname.
        line = readline(file)
        line = split(line, ':')[2]
        kname = split(line)
        readline(file)
    end

    T = readdlm(file; comments=true, comment_char='#')
    close(file)

    kline = sort(T[1:Nk, 1])
    band = Matrix{eltype(T)}(undef, Nband, Nk)
    for i in 1:Nband
        index = sortperm(T[(i-1)*Nk+1:i*Nk, 1])
        band[i, index] = T[(i-1)*Nk+1:i*Nk, 2]
    end

    if smoothband[1] ∈ ['Y', 'y']
        smoothband(band)
    end

    return band, Kline(; line = kline, name = kname, index = kindex)
end
function WAVE(file::AbstractString; orbindex="all", bandindex="all")
    error("To be continued.")
    #Read file(.dat) & file.

    datfile = open(joinpath(file * ".dat"), "r")
    readline(datfile)

    line = readline(datfile)
    wavefile = match(r"^wavefile:\s*(\S+)$", line).captures[1]

    line = readline(datfile)
    line = split(line, ':')[2]
    line = split(line, ',')
    norb = parse(Int, split(line[1], '=')[2])
    nband = parse(Int, split(line[2], '=')[2])

    line = readline(datfile)
    wavetype = match(r"^wavetype:\s*(\S+)$", line).captures[1]
    wavetype = eval(Symbol(wavetype))


    close(datfile)

    if orbindex == "all"
        orbindex = 1:norb
    end
    if bandindex == "all"
        bandindex = 1:nband
    end

    #wave
    wave = open(file, "r") do file
        return read(file)
    end

    wave = reshape(reinterpret(wavetype, wave), norb, nband)

    return wave[orbindex, bandindex]
end
"""
    smoothband(band::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}

Just as its name implies.
"""
function smoothband(band::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}

    band = deepcopy(band)
    Nk = size(band, 2)

    for i in 2:Nk
        times = 0 #A counter for avoiding dead loops.
        while true
            times += 1

            dE = band[:, i] - band[:, i-1]
            dE = dE[.!isnan.(dE)]

            (Emax, index) = findmax(abs, dE)
            if Emax < 0.2
                break
            end
            if dE[index] > 0
                band[:, i] = [NaN; band[1:end-1, i]]
            elseif dE[index] < 0
                band[:, i] = [band[2:end, i]; NaN]
            end

            if times > 100
                error("Wrong band from smoothband!")
            end
        end
    end

    I = .!any(ϵ -> isnan(ϵ), band, dims=2)[:]
    band = band[I, :]

    return band
end
function smoothband(band::Matrix{Float64}, wave::Matrix{Float64}, nk::Int)
    error("To be continued.")
    #shift some kpoints' energys to become smooth.
    Nk = size(band, 2)
    for i in 2:Nk
        times = 0 #A counter for avoiding dead loops.
        if i == nk
            (n, _) = size(wave)
            while true
                times += 1

                dE = band[:, i] - band[:, i-1]
                dE = dE[.!isnan.(dE)]

                (Emax, index) = findmax(abs, dE)
                if Emax < 0.2
                    break
                end
                if dE[index] > 0
                    band[:, i] = [NaN; band[1:end-1, i]]
                    global wave = [NaN * ones(n) wave[:, 1:end-1]]
                elseif dE[index] < 0
                    band[:, i] = [band[2:end, i]; NaN]
                    global wave = [wave[:, 2:end] NaN * ones(n)]
                end

                if times > 100
                    error("Wrong band from smoothband!")
                end
            end
        elseif i ≠ nk
            while true
                times += 1

                dE = band[:, i] - band[:, i-1]
                dE = dE[.!isnan.(dE)]

                (Emax, index) = findmax(abs, dE)
                if Emax < 0.2
                    break
                end
                if dE[index] > 0
                    band[:, i] = [NaN; band[1:end-1, i]]
                elseif dE[index] < 0
                    band[:, i] = [band[2:end, i]; NaN]
                end

                if times > 100
                    error("Wrong band from smoothband!")
                end
            end
        end
    end
    I = .!any(isnan.(band), dims=2)[:]
    band = band[I, :]
    wave = wave[:, I]
    return band, wave
end
"""
    wannierband(file::AbstractString, nk::Integer)

Read wannier90_band.dat
To be continued.
"""
function wannierband(file::AbstractString, nk::Integer)

    error("To be continued.")
    file = open(file, "r")

    kline = Vector{Float64}(undef, nk)
    E = Matrix{Float64}(undef, 0, nk)

    while !eof(file)
        dE = Vector{Float64}(undef, nk)
        for i in 1:nk
            (kline[i], dE[i]) = parse.(Float64, split(readline(file)))
        end
        E = [E; transpose(dE)]
        readline(file)
    end
    close(file)

    return E, kline
end