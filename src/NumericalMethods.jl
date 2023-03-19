module NumericalMethods

using ForwardDiff, LinearAlgebra

include("Utils.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")

using .Roots, .Min

export numderiv_one_side, numderiv_two_side, numderiv_second, gradient, hessian
export Roots.newton, Roots.secant, Roots.bisect, Roots.brent
export Min.newton, Min.brent

end # module