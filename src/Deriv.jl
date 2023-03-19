"""
    numderiv_one_side(f::Function, x::Real; δ::Real=1.0e-6)

Compute the derivative of `f` at `x` using a one-sided difference quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::Real`: The point at which to differentiate `f`.
- `δ::Real=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `deriv`: The derivative of `f` at `x`.
"""
function numderiv_one_side(f::Function, x::Real; δ::Real=1.0e-6)
    h = δ * x
    deriv = (f(x+h) - f(x))/h
    return deriv
end

"""
    numderiv_two_side(f::Function, x::Real; δ::Real=1.0e-6)

Compute the derivative of `f` at `x` using a two-sided difference quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::Real`: The point at which to differentiate `f`.
- `δ::Real=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `deriv`: The derivative of `f` at `x`.
"""
function numderiv_two_side(f::Function, x::Real; δ::Real=1.0e-6)
    h = δ * x
    deriv = (f(x+h) - f(x-h))/(2*h)
    return deriv
end


"""
    numderiv_second(f::Function, x::Real; δ::Real=1.0e-6)

Compute the second derivative of `f` at `x` using a two-sided difference
quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::Real`: The point at which to differentiate `f`.
- `δ::Real=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `deriv`: The second derivative of `f` at `x`.
"""
function numderiv_second(f::Function, x::Real; δ::Real=1.0e-6)
    first(z) = numderiv_two_side(f, z; δ=δ)
    return numderiv_two_side(first, x; δ=δ)
end

"""
    numderiv_partial(f::Function, x::AbstractVector, i::Integer; δ::Real=1.0e-6)

Compute the partial derivative of `f` at `x` with respect to the `i`th
component of `x` using a two-sided difference quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::AbstractVector`: The point at which to differentiate `f`.
- `i::Integer`: The index of the component of `x` with respect to which to
  differentiate `f`.
- `δ::Real=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `deriv`: The partial derivative of `f` at `x` with respect to the
  `i`th component of `x`.
"""
function numderiv_partial(
        f::Function, x::AbstractVector, i::Integer; δ::Real=1.0e-6
    )
    h = δ * x[i]
    h_vct = zeros(length(x)); h_vct[i] = h
    deriv = (f(x .+ h_vct) - f(x .- h_vct))/(2 * h)
    return deriv
end


"""
    gradient(f::Function, x::AbstractVector; δ::Real=1.0e-6)

Compute the gradient of `f` at `x` using a two-sided difference quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::AbstractVector`: The point at which to differentiate `f`.
- `δ::Real=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `grad::Vector`: The gradient of `f` at `x`.
"""
function gradient(f::Function, x::AbstractVector; δ::Real=1.0e-6)::Vector
    n = length(x)
    grad = Vector{eltype(x)}(undef, n)
    for i in 1:n
        grad[i] = numderiv_partial(f, x, i; δ=δ)
    end
    return grad
end

"""
    hessian(f::Function, x::AbstractVector; δ::Real=1.0e-6)

Compute the Hessian of `f` at `x` using a two-sided difference quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::AbstractVector`: The point at which to differentiate `f`.
- `δ::Real=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `hess::Matrix`: The Hessian of `f` at `x`.
"""
function hessian(f::Function, x::AbstractVector; δ::Real=1e-6)::Matrix
    n = length(x)
    hess = Matrix{eltype(x)}(undef, n, n)
    for i ∈ 1:n, j ∈ 1:n
        partial_i(z) = numderiv_partial(f, z, i; δ=δ)
        hess[i, j] = numderiv_partial(partial_i, x, j; δ=δ)
    end
    return hess
end