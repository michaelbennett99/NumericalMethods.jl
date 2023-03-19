module NumericalMethods

using ForwardDiff, LinearAlgebra
using .Roots, .Min

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")

end # module