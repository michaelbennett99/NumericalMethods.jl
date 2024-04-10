@testset "Roots" begin
    f(x) = exp(x-1.5) - (x-1.5)^2 - 5
    f_root = 3.856

    @testset "bisect" begin
        for low in -10:1:3
            for high in 5:1:20
                x, fx, _ = Roots.bisect(f, low, high)
                @test x ≈ f_root
                @test fx ≈ 0.0
            end
        end
        @test_throws ArgumentError Roots.bisect(f, -10., 3.)
        @test_throws ArgumentError Roots.bisect(f, 5., 10.)
    end

    @testset "secant" begin
        for x_0 in -10.:1.:10.
            for x_1 in -10.:1.:10.
                if (x_0 != x_1)
                    x, fx, _ = Roots.secant(f, x_0, x_1)
                    @test x ≈ f_root
                    @test fx ≈ 0.0
                end
            end
        end
        @test_throws ArgumentError Roots.secant(f, 0., 0.)
    end

    @testset "newton" begin
        for x_0 in -10.:1.:10.
            x, fx, _ = Roots.newton(f, x_0)
            @test x ≈ f_root
            @test fx ≈ 0.0
        end
    end

    @testset "brent" begin
        for low in -10.:1.:3.
            for high in 5.:1.:20.
                x, fx, _ = Roots.brent(f, low, high)
                @test x ≈ f_root
                @test fx ≈ 0.0
            end
        end
        @test_throws ArgumentError Roots.brent(f, -10., 3.)
        @test_throws ArgumentError Roots.brent(f, 5., 10.)
    end
end
