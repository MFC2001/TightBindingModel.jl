
function Insulator_distribution(Ne::Integer, N::Integer)::Function
	T = [ones(Int, Ne); zeros(Int, N - Ne)]
	f(x) = T
	return f
end

function FDdistribution(μ::Real, T::Real)::Function
	#energy unit is eV, so use kB/e.
	kB = 8.6173332621e-5

	kBT = kB * T

	if T == 0
		f0(ϵ) = ϵ > μ ? 0 : 1
		return f0
	elseif T > 0
		f(ϵ) = 1 / (exp((ϵ - μ) / kBT) + 1)
		return f
	else
		error("Tempreture should be not less than 0. ")
	end
end