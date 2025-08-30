export atom_orbindex, atom2orbindex
"""
    atom_orbindex(orbital::ORBITAL, poscar::POSCAR; deps::Real=0.5, spinindex::AbstractVector{<:Pair}=Pair[])
"""
function atom_orbindex(orbital::ORBITAL, poscar::POSCAR; deps::Real=0.5, spinindex::AbstractVector{<:Pair}=Pair[])

    error("To be continued.")
    #Classify all orbitals of unitcell.
    n = poscar.num
    atomlocat = poscar.location
    deps = deps^2

    if isempty(spinindex)
        orblocat = orbital.location

        orbital_index = Vector{Vector{Int}}(undef, n)
        Tindex = collect(1:orbital.num)
        for i in 1:n
            D = sum(abs2, orblocat .- transpose(atomlocat[i, :]), dims=2)[:]
            I = D .< deps
            orbital_index[i] = Tindex[I]
            I = .!I
            orblocat = orblocat[I, :]
            Tindex = Tindex[I]
        end

        if sum(length.(orbital_index)) ≠ orbital.num
            error("The ORBITAL or POSCAR may be incompatible.")
        end

    else
        orbindex = orbital.index
        orbital_spin_index = Vector{Pair}(undef, 0)
        for si in spinindex
            T = Vector{Int}(undef, 0)
            for i in eachindex(orbindex)
                if orbindex[i] ∈ si.second
                    push!(T, i)
                end
            end
            push!(orbital_spin_index, si.first => T)
        end

        orbital_index = Vector{Pair}(undef, 0)
        for si in orbital_spin_index
            Tindex = si.second
            orblocat = orbital.location[Tindex, :]
            Torbital_index = Vector{Vector{Int}}(undef, n)
            for i in 1:n
                D = sum(abs2, orblocat .- transpose(atomlocat[i, :]), dims=2)[:]
                I = D .< deps
                Torbital_index[i] = Tindex[I]
                I = .!I
                orblocat = orblocat[I, :]
                Tindex = Tindex[I]
            end

            if sum(length.(Torbital_index)) ≠ length(si.second)
                error("The ORBITAL or POSCAR may be incompatible.")
            end

            push!(orbital_index, si.first => Torbital_index)
        end

    end

    return orbital_index
end
"""
    atom2orbindex(atom_index::AbstractVector{<:Pair}, atom_orbital_index::AbstractVector)

"""
function atom2orbindex(atom_index::AbstractVector{<:Pair}, atom_orbital_index::AbstractVector)
    orbital_index = similar(atom_index, Pair{String,Vector{Int}})
    for i in eachindex(atom_index)
        orbindex = Vector{Int}(undef, 0)
        for j in atom_index[i].second
            append!(orbindex, atom_orbital_index[j])
        end
        orbital_index[i] = atom_index[i].first => orbindex
    end
    return orbital_index
end