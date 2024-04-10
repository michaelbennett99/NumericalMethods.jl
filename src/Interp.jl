module Interp

using StaticArrays, LinearAlgebra

export AbstractInterpolator, linear_interpolation, cubic_interpolation

abstract type AbstractInterpolator end

"""
A callable linear interpolation object.
"""
struct Linear{N, T} <: AbstractInterpolator
    x::SVector{N, T}
    y::SVector{N, T}
end

"""
    linear_interpolation(x, y)

Construct a linear interpolation function from the given points.

# Arguments

- `x::Vector{<:Real}`: The x-coordinates of the points.
- `y::Vector{<:Real}`: The y-coordinates of the points.

# Returns

- `f::Linear`: A callable linear interpolation object.
"""
function linear_interpolation(x::Vector{<:Real}, y::Vector{<:Real})::Linear
    if ! issorted(x)
        error("x must be sorted")
    end
    return Linear{length(x), eltype(x)}(x, y)
end

"""
    (f::Linear)(x)

Evaluate the linear interpolation function at `x`.

# Arguments

- `f::Linear`: A callable linear interpolation object.
- `x::Real`: The point at which to evaluate the interpolation function.

# Returns

- `y::Real`: The value of the interpolation function at `x`.
"""
function (f::Linear)(x)
    if x < f.x[1] || x > f.x[end]
        error("x is out of range")
    elseif x ∈ f.x
        return f.y[findfirst(f.x .== x)]
    else
        lub = findfirst(f.x .> x)
        glb = lub - 1
        @views grad = (f.y[lub] - f.y[glb]) / (f.x[lub] - f.x[glb])
        return f.y[glb] + grad * (x - f.x[glb])
    end
end


"""
A callable cubic spline interpolation object.
"""
struct Spline{N, T} <: AbstractInterpolator
    x::SVector{N, T}
    y::SVector{N, T}
    k::SVector{N, T}
end

"""
    cubic_interpolation(x, y)

Construct a cubic interpolation function from the given points.

# Arguments

- `x::Vector{<:Real}`: The x-coordinates of the points.
- `y::Vector{<:Real}`: The y-coordinates of the points.

# Returns

- `f::Spline`: A callable cubic interpolation object.
"""
function cubic_interpolation(
        x::Vector{<:Real}, y::Vector{<:Real}
    )::Spline{length(x), Float64}
    if ! issorted(x)
        error("x must be sorted")
    end
    n = length(x)
    dl = zeros(n-1)
    d = zeros(n)
    du = zeros(n-1)
    r = zeros(n)

    dl[1] = 1/(x[2] - x[1])
    du[1] = 1/(x[2] - x[1])
    d[1] = 2/(x[2] - x[1])
    d[n] = 2/(x[n] - x[n-1])
    r[1] = 3 * (y[2] - y[1]) / (x[2] - x[1])^2
    r[n] = 3 * (y[n] - y[n-1]) / (x[n] - x[n-1])^2

    for i ∈ 2:n-1
        dx1 = x[i] - x[i-1]
        dx2 = x[i+1] - x[i]
        dy1 = y[i] - y[i-1]
        dy2 = y[i+1] - y[i]
        d[i] = 2 * (1/dx1 + 1/dx2)
        dl[i] = 1/dx2
        du[i] = 1/dx2
        r[i] = 3 * (dy1/dx1^2 + dy2/dx2^2)
    end

    vdp = Tridiagonal(dl, d, du) \ r
    return Spline{n, Float64}(x, y, vdp)
end

"""
    (f::Spline)(x)

Evaluate the cubic interpolation function at `x`.

# Arguments

- `f::Spline`: A callable cubic interpolation object.
- `x::Real`: The point at which to evaluate the interpolation function.

# Returns

- `y::Real`: The value of the interpolation function at `x`.
"""
function (f::Spline)(x::Real)
    if x < f.x[1] || x > f.x[end]
        error("x is out of range")
    elseif x ∈ f.x
        return f.y[findfirst(f.x .== x)]
    else
        lub = findfirst(f.x .> x)
        if lub === nothing
            println(x)
            println(f.x)
            error("x is out of range")
        end
        glb = lub - 1
        @views t = (x - f.x[glb]) / (f.x[lub] - f.x[glb])
        @views a = f.k[glb] * (f.x[lub] - f.x[glb]) - (f.y[lub] - f.y[glb])
        @views b = -f.k[lub] * (f.x[lub] - f.x[glb]) + (f.y[lub] - f.y[glb])
        @views y_hat = (((1 - t) * f.y[glb]) + (t * f.y[lub])
            + (t * (1 - t) * (a * (1 - t) + b * t)))
        return y_hat
    end
end

end # module
