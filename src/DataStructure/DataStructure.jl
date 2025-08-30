

const Mat3{T} = SMatrix{3, 3, T, 9}
const Vec3{T} = SVector{3, T}

include("./CrystallographyCore/CrystallographyCore.jl")
include("./Hopping.jl")
include("./HR.jl")
include("./ORBITAL.jl")
include("./AbstractReciprocalHoppings/AbstractReciprocalHoppings.jl")
include("./AbstractInterAction/AbstractInterAction.jl")

include("./WannierBasis.jl")
include("./SymOp.jl")
include("./BrillouinZone.jl")

include("./vckArray.jl")
include("AbstractModel/AbstractModel.jl")
