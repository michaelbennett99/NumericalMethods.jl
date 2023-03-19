module NumericalMethods

using ForwardDiff, LinearAlgebra
using .Roots, .Min

export numderiv_one_side, numderiv_two_side, numderiv_second, gradient, hessian
export bisect, secant, brent, newton

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")

end # module