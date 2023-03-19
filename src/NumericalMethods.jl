module NumericalMethods

using ForwardDiff, LinearAlgebra

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")

using .Deriv, .Roots, .Min

export Deriv, Roots, Min

end # module