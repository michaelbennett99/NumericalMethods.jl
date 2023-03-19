module NumericalMethods

using ForwardDiff, LinearAlgebra

export numderiv_one_side, numderiv_two_side, numderiv_second, gradient, hessian
export Roots.bisect, Roots.secant, Roots.brent, Roots.newton
export Min.brent, Min.newton

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")

end # module