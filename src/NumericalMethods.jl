module NumericalMethods

using ForwardDiff, LinearAlgebra

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")
include("Interp.jl")

using .Deriv, .Roots, .Min, .Interp

export Deriv, Roots, Min, Interp

end # module