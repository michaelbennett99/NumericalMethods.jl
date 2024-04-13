module NumericalMethods

using ForwardDiff, LinearAlgebra

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")
include("Interp.jl")

using .Deriv
using .Roots
using .Min
using .Interp

export Deriv
export Roots
export Min
export Interp

end # module
