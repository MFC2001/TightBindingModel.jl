struct AtomPath
	unitpath::SVector{3, Int} #path of unit lattice vector
	hrindex::Vector{Int}     #the index of HR
end
struct AtomHR
	parindex::Matrix{Tuple{Int, Int}}  #participant atom index, index[1] to index[2], [n, n]
	parname::Matrix{String} #participator atom name, "name[1]->name[2]", [n, n]
	atompath::Matrix{Vector{AtomPath}} #[n, n][N],[ [Atompath] ]
	Natompath::Matrix{Int} #[n, n],[ length([Atompath]) ]
	maxlattdist::Float64 #the max distance of unitcell path.
	hrorbindex::Matrix{Int} #orbital_index corresponding to hr.path, reindex of each atom.
	orbindex::Vector{Vector{Int}} #orbital index corresponding to atom.
	num_orb::Vector{Int} #number of each atom's orbitals.
end
