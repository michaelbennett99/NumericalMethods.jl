module NumericalMethods

using ForwardDiff, LinearAlgebra

export numderiv_one_side, numderiv_two_side, numderiv_second
export gradient, hessian
export bisect, secant, func_iter, newton, brent

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")

end # module