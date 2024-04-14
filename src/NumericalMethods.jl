module NumericalMethods

using ForwardDiff, LinearAlgebra

include("Errors.jl")
include("Deriv.jl")
include("Roots.jl")
include("Min.jl")
include("Interp.jl")

export Deriv
export Roots
export Min
export Interp

end # module
