export blockdiagm, blockdiagm!
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


    T = zeros(eltype(args[1].second[1]), a, b)

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