using NumericalMethods
using Test

@testset "Differentiation" begin
    x = collect(range(0, 2π, step=π/4))
    for x_i in x
        @test numderiv_one_side(sin, x_i) ≈ cos(x_i)
        @test numderiv_two_side(sin, x_i) ≈ cos(x_i)
        @test numderiv_second(sin, x_i) ≈ -sin(x_i)
    end
    f(x, y) = x^2 + y^2
    f_x(x, y) = 2x
    f_y(x, y) = 2y
    f_xx(x, y) = 2
    J(x, y) = [f_x(x, y) f_y(x, y)]
    H(x, y) = [f_xx(x, y) 0; 0 f_xx(x, y)]
    for x_i in x
        for y_i in x
            @test gradient(f, [x_i, y_i]) ≈ J(x_i, y_i)
            @test hessian(f, [x_i, y_i]) ≈ H(x_i, y_i)
        end
    end
    # Write your tests here.
end
