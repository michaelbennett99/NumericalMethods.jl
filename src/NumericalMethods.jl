module NumericalMethods

using ForwardDiff, LinearAlgebra

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")

using .Roots, .Min

print(typeof(Roots.newton))

end # module