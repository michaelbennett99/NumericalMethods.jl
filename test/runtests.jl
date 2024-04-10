using Test
using NumericalMethods

≈(x, y) = isapprox(x, y; atol=1e-4, rtol=1e-4)

include("Deriv.jl")
include("Interp.jl")
include("Min.jl")
include("Roots.jl")
