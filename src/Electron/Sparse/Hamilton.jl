function SparseHamilton(k::AbstractVector{<:Real}, hr::HR, lattvec::AbstractMatrix{<:Real})
	norb = hr.norb
	kl = transpose(lattvec) * k

	H = spzeros(ComplexF64, norb, norb)

	for j in 1:norb
		for i in 1:j #Access by column.
			if hr.Nindex[i, j] > 0
				H[i, j] = sum(hr.index[i, j]) do index
					return hr.value[index] * cis(kl ⋅ hr.path[index])
				end
			end
		end
	end

	#The Hermitian is just a shell.
	return Hermitian(H, :U)
end
