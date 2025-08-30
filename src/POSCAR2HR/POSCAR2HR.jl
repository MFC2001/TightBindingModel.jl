module POSCAR2HRs

using LinearAlgebra
using StaticArrays
import ..LatticeModel: Lattice, Cell, ReducedCoordinates, CartesianCoordinates, ORBITAL, HR

include("./AtomHR.jl")
include("./SuperCell.jl")
include("./check.jl")
include("./findunitcell.jl")
include("./hrsplit.jl")
include("./others.jl")
include("./postprocess.jl")

include("./atompath_sum.jl")
include("./ExpandHR.jl")

include("./ContractHR.jl")

end

import .POSCAR2HRs: POSCAR2HRs, ExpandHR, ContractHR
export ExpandHR, ContractHR

include("./shell.jl")
