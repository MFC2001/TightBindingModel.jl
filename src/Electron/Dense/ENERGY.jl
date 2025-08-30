export ENERGY
function ENERGY(k::AbstractVector{<:Real}, args...; vector = false)

	H = Hamilton(k, args...)

	if vector
		return eigen!(H)
	else
		return eigvals!(H)
	end
end

function ENERGY(hr::HR; vector = false)

	H = Hamilton(hr)

	if vector
		return eigen!(H)
	else
		return eigvals!(H)
	end
end

# function ENERGY(hr::HR, Δ::Number; vector=false)

#     H = HamiltonSC(hr, Δ)

#     if vector
#         return eigen!(H)
#     else
#         return eigvals!(H)
#     end
# end
