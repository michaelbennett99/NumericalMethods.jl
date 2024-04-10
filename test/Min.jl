@testset "Min" begin
    f(x) = (x - 2)^2
    f_argmin = 2.0
    f_min = 0.0

    g(x) = (x[1] - 3)^2 + (x[2] - 4)^2
    g_argmin = [3.0, 4.0]
    g_min = 0.0

    @testset "Brent's Method" begin
        test_params = [
            (0.0, 0.5, 5.0),
            (-3.0, 2.0, 10.0),
            (1.0, 1.5, 2.5),
            (-10.0, 0.0, 10.0)
        ]
        for (a, b, c) in test_params
            brent_min = Min.brent(f, a, b, c)
            @test brent_min[1] ≈ f_argmin
            @test brent_min[2] ≈ f_min
        end
    end

    @testset "Newton's Method Univariate" begin
        for x_0 ∈ [-5.0, 0.0, 1.0, 5.0, 10.0]
            for δ ∈ [1e-4, 1e-6, 1e-8]
                newton_min = Min.newton(f, x_0, δ)
                @test newton_min[1] ≈ 2.0
                @test newton_min[2] ≈ 0.0
            end

            newton_auto = Min.newton(f, x_0)
            @test newton_auto[1] ≈ f_argmin
            @test newton_auto[2] ≈ f_min
        end
    end

    @testset "Newton's Method Multivariate" begin
        for x_0 ∈ [[-5.0, -5.0], [0.0, 0.0], [1.0, 1.0], [5.0, 5.0], [10.0, 10.0]]
            for δ ∈ [1e-4, 1e-6, 1e-8]
                newton_min = Min.newton(g, x_0, δ)
                @test newton_min[1] ≈ g_argmin
                @test newton_min[2] ≈ g_min
            end

            newton_auto = Min.newton(g, x_0)
            @test newton_auto[1] ≈ g_argmin
            @test newton_auto[2] ≈ g_min
        end
    end
end
