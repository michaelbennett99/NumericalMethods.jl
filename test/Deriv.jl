using NumericalMethods
using TestItems

@testitem "Univariate" begin
    ≈(x, y) = isapprox(x, y; atol=1e-4, rtol=1e-4)
    x = collect(range(0, 2π, step=π/4))
    for x_i in x
        @test Deriv.differentiate(sin, x_i; two_side=false) ≈ cos(x_i)
        @test Deriv.differentiate(sin, x_i) ≈ cos(x_i)
        @test Deriv.twice_differentiate(sin, x_i) ≈ -sin(x_i)
    end
end

@testitem "Multivariate" begin
    ≈(x, y) = isapprox(x, y; atol=1e-4, rtol=1e-4)
    x = collect(range(0, 2π, step=π/4))
    f(x) = x[1]^2 + x[2]^2
    f_1(x) = 2x[1]
    f_2(x) = 2x[2]
    f_zz(x) = 2
    J(x) = [f_1(x), f_2(x)]
    H(x) = [f_zz(x) 0; 0 f_zz(x)]
    for x_i in x
        for y_i in x
            z_i = [x_i, y_i]
            @test Deriv.gradient(f, z_i) ≈ J(z_i)
            @test Deriv.hessian(f, z_i) ≈ H(z_i)
        end
    end
end