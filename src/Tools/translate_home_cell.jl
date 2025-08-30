export translate_home_cell
"""
    translate_home_cell(lattvec::AbstractMatrix{<:Real}, orbital::ORBITAL, hr::HR, translatepair::Pair{<:Integer,<:AbstractVector{<:Integer}}...)

    translatepair is "orbital index" => "orbital's location in positive lattice", orbital's location = orbital.location["orbital index"].
"""
function translate_home_cell(lattvec::AbstractMatrix{<:Real}, orbital::ORBITAL, hr::HR, translatepair::Pair{<:Integer,<:AbstractVector{<:Integer}}...)

    neworbital = translate_home_cell_orbital(lattvec, orbital, translatepair...)
    newhr = translate_home_cell_hr(hr, translatepair...)

    return neworbital, newhr
end
function translate_home_cell_hr(hr::HR, translatepair::Pair{<:Integer,<:AbstractVector{<:Integer}}...)

    newhr = deepcopy(hr)

    for p in translatepair

        for i in setdiff(1:newhr.norb, p.first)
            if newhr.Nindex[i, p.first] > 0
                for index in newhr.index[i, p.first]
                    newhr.path[index] += p.second
                end
            end
        end
        for i in setdiff(1:newhr.norb, p.first)
            if newhr.Nindex[p.first, i] > 0
                for index in newhr.index[p.first, i]
                    newhr.path[index] -= p.second
                end
            end
        end

    end

    return newhr
end
function translate_home_cell_orbital(lattvec::AbstractMatrix{<:Real}, orbital::ORBITAL, translatepair::Pair{<:Integer,<:AbstractVector{<:Integer}}...)

    neworbital = deepcopy(orbital)

    for p in translatepair
        neworbital.location[p.first] -= lattvec * p.second
    end

    return neworbital
end