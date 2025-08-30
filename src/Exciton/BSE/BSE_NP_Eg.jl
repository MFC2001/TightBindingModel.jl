
function BSE_NP_Eg(BSEband, vckindex, bandk, bandkq, scissor)

	ittr = eachindex(vckindex)

	Eg = function (n)
        A = BSEband.vectors[:, n]
		return sum(ittr) do i
			(v, c, k) = vckindex[i]
			abs2(A[i]) * (bandkq[k].values[c] - bandk[k].values[v] + scissor)
		end
	end

	return Eg
end
