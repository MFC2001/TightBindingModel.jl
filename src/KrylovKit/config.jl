
using Main.KrylovKit
using Main.KrylovDefaults

export eigsolveconfigure!

const eigsolveconfig = Dict{Symbol, Any}(
	:howmany => 1,
	:which => EigSorter(abs; rev = false),
	:verbosity => 0,
	:tol => KrylovDefaults.tol,
	:krylovdim => KrylovDefaults.krylovdim,
	:maxiter => 300,
	:orth => KrylovDefaults.orth,
	:issymmetric => false,
	:ishermitian => false,
	:eager => false,
)

function eigsolveconfigure!(; kwards...)
	for (k, v) in kwards
		eigsolveconfig[k] = v
		if k == :howmany
			eigsolveconfig[:krylovdim] = max(2 * v + 10, 30)
		end
	end
end
