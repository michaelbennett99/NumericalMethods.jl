f(x) = 2 * x^2 + 3 * x + 1

@testset "Linear Interpolation" begin
    x = [5, 6]
    y = f.(x)
    linear_interpolation = Interp.linear_interpolation(x, y)

    @test linear_interpolation(5) == y[1]
    @test linear_interpolation(5.5) ≈ 0.5 * y[1] + 0.5 * y[2]
    @test_throws ErrorException linear_interpolation(0)

    wrong_x = [6, 5]
    wrong_y = f.(wrong_x)
    @test_throws ErrorException Interp.linear_interpolation(wrong_x, wrong_y)
end

@testset "Cubic Interpolation" begin
    x = [5, 6, 7, 8, 9]
    y = f.(x)
    cubic_interpolation = Interp.cubic_interpolation(x, y)

    @test cubic_interpolation(5) == y[1]
    @test isapprox(cubic_interpolation(5.5), f(5.5); atol=1e-2, rtol=1e-2)
    @test cubic_interpolation(6) == y[2]
    @test isapprox(cubic_interpolation(6.5), f(6.5); atol=1e-2, rtol=1e-2)
    @test cubic_interpolation(7) == y[3]
    @test_throws ErrorException cubic_interpolation(0)

    wrong_x = [6, 5, 7]
    wrong_y = f.(wrong_x)
    @test_throws ErrorException Interp.cubic_interpolation(wrong_x, wrong_y)
end
