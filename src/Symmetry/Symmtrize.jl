export Symmtrize
"""

"""
function Symmtrize(TB::TightBindModel{HT}, symops::Vector{SymOp}; heps = 1e-2, reps = 1e-3) where {HT <: Real}

	# TODO without spin presently.


	I = CartesianIndices(TB.hop)
	symhop = [Hopping{HT}[] for _ in I]
	hopmap = [Function[] for _ in I]
	irredhop = Hopping{Int}[]
	symvalue = HT[]

	hrindex = Matrix{Vector{Int}}(undef, size(TB.hop)...)
	allhop = Hopping{HT}[]
	for II in I
		hrindex[II] = collect(length(allhop) .+ eachindex(TB.hop[II]))
		append!(allhop, TB.hop[II])
	end


	index = 0
	while true
		hop = _get_first_hop(allhop, hrindex)
		if isnothing(hop)
			break
		end
		index += 1
		push!(irredhop, Hopping{Int}(hop.i, hop.j, hop.R, 0))
		push!(symvalue, hop.t)

		for symop in symops
			equivlent_hop = symop_hop(symop, hop, TB.orb_location; symtol = reps)

			i, j = equivlent_hop.i, equivlent_hop.j
			I = findfirst(x -> x.R == equivlent_hop.R && isapprox(x.t, equivlent_hop.t; atol = heps), allhop[hrindex[i, j]])
			if !isnothing(I)
				deleteat!(hrindex[i, j], I)
				push!(symhop[i, j], equivlent_hop)
				let i = i, j = j, R = equivlent_hop.R, index = index
					push!(hopmap[i, j], t -> Hopping{T₁}(i, j, R, t[index]))
					# t[index])
				end
				irredhop[end].t += 1
			end
			# Hermitian.
			I = findfirst(x -> x.R == -equivlent_hop.R && isapprox(x.t, equivlent_hop.t; atol = heps), allhop[hrindex[j, i]])
			if !isnothing(I)
				deleteat!(hrindex[j, i], I)
				push!(symhop[j, i], conj(equivlent_hop))
				let i = i, j = j, R = equivlent_hop.R, index = index
					push!(hopmap[j, i], t -> Hopping{T₁}(i, j, R, conj(t[index])))
				end
				irredhop[end].t += 1
			end
		end
	end

	return SymTightBindModel{HT}(TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location,
		symhop, TB.Nhop, TB.period, symops, hopmap, irredhop, symvalue)
end


function symop_hop(symop, hop::Hopping{T}, orblocation; symtol = 1e-6) where {T}
	startpoint = symop.W * orblocation[hop.i] + symop.w
	endpoint = symop.W * (orblocation[hop.j] + hop.R) + symop.w

	startorbital = findfirst(x -> _is_approx_integer(startpoint - x; atol = symtol), orblocation)
	endorbital = findfirst(x -> _is_approx_integer(endpoint - x; atol = symtol), orblocation)

	path = round.((endpoint - orblocation[endorbital]) - (startpoint - orblocation[startorbital]))

	return Hopping{T}(startorbital, endorbital, path, hop.t)
end


function _get_first_hop(allhop, hrindex)
	for index in hrindex
		if !isempty(index)
			return allhop[index[1]]
		end
	end
	return nothing
end
