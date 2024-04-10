using Test
using NumericalMethods

≈(x, y) = isapprox(x, y; atol=1e-4, rtol=1e-4)

@testset "Deriv" begin
    include("Deriv.jl")
end
