module NumericalMethods

using ForwardDiff, LinearAlgebra

export numderiv_one_side, numderiv_two_side
export bisect, secant, func_iter, newton, brent

"""
    numderiv_one_side(f::Function, x::Real; δ::Real=1.0e-6)

Compute the derivative of `f` at `x` using a one-sided difference quotient.

# Arguments

- `f::Function`: The function to differentiate.
- `x::T`: The point at which to differentiate `f`.
- `δ::T=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `deriv::T`: The derivative of `f` at `x`.
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
- `x::T`: The point at which to differentiate `f`.
- `δ::T=1.0e-6`: The step size to use in the difference quotient.

# Returns

- `deriv::T`: The derivative of `f` at `x`.
"""
function numderiv_two_side(f::Function, x::Real; δ::Real=1.0e-6)
    h = δ * x
    deriv = (f(x+h) - f(x-h))/(2*h)
    return deriv
end

function numderiv_partial(
        f::Function, x::AbstractVector, i::Integer; δ::Real=1.0e-6
    )
    h = δ * x
    h_vct = zeros(len(x)); h_vct[i] = h
    deriv = (f(x+h_vct) - f(x-h_vct))/(2 * h)
    return deriv
end

function gradient(f::Function, x::AbstractVector; δ::Real=1.0e-6)
    n = length(x)
    grad = Matrix{Float64}(undef, n, 1)
    for i in 1:n
        grad[i, 1] = numderiv_partial(f, x, i; δ=δ)
    end
    return grad
end

function hessian(f::Function, x::AbstractVector; δ::Real=1e-6)
    n = length(x)
    hess = Matrix{Float64}(undef, n, n)
    for i ∈ 1:n, j ∈ 1:n
        partial_i(z) = numderiv_partial(f, z, i; δ=δ)
        hess[i, j] = numderiv_partial(partial_i, x, j; δ=δ)
    end
    return hess
end

"""
    bisect(f, xlow, xhigh; tol=1.0e-6, maxiter=30, verbose=false, kwargs...)

Find the root of `f` in the interval `[xlow, xhigh]` using the bisection method.

Such a method will not work in multiple dimensions.

# Arguments

- `f::Function`: The function to find the root of.
- `xlow::Real`: The lower bound of the interval.
- `xhigh::Real`: The upper bound of the interval.
- `tol::Real=1.0e-6`: The tolerance for the root.
- `maxiter::Integer=30`: The maximum number of iterations to perform.
- `verbose::Bool=false`: Whether to print the iteration history.

# Returns

- `xroot::Real`: The root of `f`.
- `fxroot::Real`: The value of `f` at the root.
- `xlow::Real`: The lower bound of the root interval.
- `xhigh::Real`: The upper bound of the root interval.
- `niter::Integer`: The number of iterations performed.
"""
function bisect(
        f::Function, xlow::Real, xhigh::Real;
        tol::Real=1.0e-6, maxiter::Integer=30, verbose::Bool=false, kwargs...
    )
    fxlow = f(xlow; kwargs...)
    fxhigh = f(xhigh; kwargs...)
    
    if (fxlow*fxhigh > 0) 
       throw(ArgumentError("Root not bracketed."))
    elseif (xlow > xhigh) 
       throw(ArgumentError("Brackets out of order.")) 
    end
    
    diff = xhigh - xlow
    downwardsloping = (fxlow > 0)
    
    niter = 0
    while niter < maxiter
        niter += 1
        xcur = (xlow + xhigh) / 2
        fxcur = f(xcur; kwargs...)
        if downwardsloping
            if (fxcur > 0) 
                xlow = xcur
            else
                xhigh = xcur
            end
        else
            if (fxcur > 0) 
                xhigh = xcur
            else 
                xlow = xcur
            end
        end               
        diff = abs(xhigh-xlow)
        
        if verbose
            writeio(
                stdout, (4, (5, 15.8)),
                niter, xlow, xhigh, xcur, fxcur, diff, callwait=false
            )
        end        
        
        if (diff < tol) 
            break
        end
    end
 
    if (diff > tol)
        throw(
            ArgumentError(
                "Did not converge: xlow = $xlow, xhigh = $xhigh, diff = $diff"
            )
        )           
    end
 
    xroot = (xlow + xhigh)/2
    fxroot = f(xroot; kwargs...)
    
    return xroot, fxroot, xlow, xhigh, niter
end

"""
    secant(f, x_0, x_1; tol=1e-6, max_iter=1000, kwargs...)

Use the secant method to find the root of a function f.

# Arguments

- `f`: The function to find the root of.
- `x_0`: The first guess.
- `x_1`: The second guess.
- `tol=1e-6`: The tolerance for the root.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
- `iter`: The number of iterations performed.
"""
function secant(
        f::Function, x_0::Real, x_1::Real;
        tol::Real=1e-6, max_iter::Integer=1000, kwargs...
    )
    iter = 0
    diff = tol + 1
    while diff > tol
        iter += 1

        sec = (x_1 - x_0) / (f(x_1; kwargs...) - f(x_0; kwargs...))
        x_2 = x_1 - f(x_1; kwargs...) * sec

        diff = abs(x_2 - x_1)

        x_0, x_1 = x_1, x_2

        if iter == max_iter
            break
        end
    end

    if (diff > tol)
        throw(ArgumentError("Did not converge: diff = $diff"))           
    end

    xroot = x_1
    fxroot = f(xroot; kwargs...)

    return xroot, fxroot, iter
end

"""
    func_iter(f, x_0; tol=1e-6, max_iter=1000, kwargs...)

Use the functional iteration method to find the root of a function f. This
method works for functions from R^n to R^n, where n ≥ 1.

# Arguments

- `f::Function`: The function to find the root of.
- `x_0`: The first guess.
- `λ=0`: The relaxation parameter.
- `tol=1e-6`: The tolerance for the root.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
"""
function func_iter(
        f::Function, x_0::Union{AbstractVector{<:Real}, Real};
        λ::Real=0, tol=1e-6, max_iter::Integer=1000, kwargs...
    )
    g(x) = f(x; kwargs...) .+ x
    
    iter = 0
    diff = tol + 1
    while diff > tol
        iter += 1

        x_1 = (1 - λ) .* g(x_0) .+ λ .* x_0

        diff = norm(x_1 - x_0)
        
        x_0 = x_1

        if iter == max_iter
            break
        end
    end

    if (diff > tol)
        throw(ArgumentError("Did not converge: diff = $diff"))           
    end

    xroot = x_0
    fxroot = f(xroot; kwargs...)

    return xroot, fxroot, iter
end


"""
    newton(f::Function, x_0; tol=1e-6, max_iter=1000, kwargs...)

Use Newton's method to find the root of a function f.

# Arguments

- `f::Function`: The function to find the root of.
- `x_0`: The first guess.
- `tol=1e-6`: The tolerance for the root.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
- `iter`: The number of iterations performed.
"""
function newton(
        func::Function, x_0::Float64;
        tol::Real=1e-6, max_iter::Integer=1000, kwargs...
    )
    f(x) = func(x; kwargs...)
    df(x) = ForwardDiff.derivative(f, x) # Using automatic differentiation
    # provided by the ForwardDiff package . See juliadiff.org for justification
    # of use of automatic differentiation
    iter = 0
    diff = tol + 1
    while diff > tol
        iter += 1

        x_1 = x_0 - f(x_0) / df(x_0)

        diff = abs(x_1 - x_0)

        x_0 = x_1

        if iter == max_iter
            break
        end
    end

    if (diff > tol)
        throw(ArgumentError("Did not converge: diff = $diff"))           
    end

    xroot = x_0
    fxroot = f(xroot)

    return xroot, fxroot, iter
end

"""
    newton(
        f::Function, x_0::AbstractVector{<:Real};
        tol::Real=1e-6, max_iter::Integer=1000, kwargs...
    )

Use Newton's method to find the root of a multivariate function f.

# Arguments

- `f::Function`: The function to find the root of.
- `x_0`: The first guess.
- `tol=1e-6`: The tolerance for the root.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
- `iter`: The number of iterations performed.
"""
function newton(
        func::Function, x_0::Union{AbstractVector{<:Real}, Real};
        tol::Real=1e-6, max_iter::Integer=1000, kwargs...
    )
    f(x) = func(x; kwargs...)
    J(x) = ForwardDiff.jacobian(f, x) # Using automatic differentiation
    # See juliadiff.org for justifications of using automatic differentiation
    iter = 0
    diff = tol + 1
    while diff > tol
        iter += 1

        x_1 = x_0 .- inv(J(x_0)) * f(x_0)

        diff = norm(x_1 - x_0)

        x_0 = x_1

        if iter == max_iter
            break
        end
    end

    if (diff > tol)
        throw(ArgumentError("Did not converge: diff = $diff"))           
    end

    xroot = x_0
    fxroot = f(xroot)

    return xroot, fxroot, iter
end

"""
    brent(f::Function, x_0, x_1; tol=1e-6, max_iter=1000, kwargs...)

Use Brent's method to find the root of a function f.

# Arguments

- `f::Function`: The function to find the root of.
- `x_0`: The first guess.
- `x_1`: The second guess.
- `rtol=1e-6`: The tolerance for the root.
- `ftol=1e-6`: The tolerance for the function value.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
"""
function brent(
        func::Function, x_0::Real, x_1::Real;
        rtol::Real=1e-6, ftol::Real=1e-6, max_iter::Integer=1000, kwargs...
    )
    f(x) = func(x; kwargs...)

    global d, e, iter

    toler = ftol
    tol1 = rtol

    a = x_0
    b = x_1
    fa = f(a)
    fb = f(b)

    if (fa*fb > 0) 
        throw(ArgumentError("Root not bracketed."))
    elseif (a > b) 
        throw(ArgumentError("Brackets out of order.")) 
    end

    c = b
    fc = fb

    for i ∈ 1:max_iter
        if (fb*fc > 0)
            c = a
            fc = fa
            d = b - a
            e = d
        end

        if (abs(fc) < abs(fb))
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        
        xm = 0.5 * (c - b)
        if (abs(xm) <= tol1 || fb == 0)
            iter = i
            break
        end

        if (abs(e) >= tol1 && abs(fa) > abs(fb))
            s = fb/fa
            if (abs(c - a) ≤ toler)
                p = 2 * xm * s
                q = 1 - s
            else
                q = fa/fc
                r = fb/fc
                p = s * (2 * xm * q * (q-r) - (b-a) * (r-1))
                q = (q-1) * (r-1) * (s-1)
            end

            if (p > 0)
                q = -q
            else
                p = -p
            end

            term1 = 2 * p
            term2 = min(3 * xm * q - abs(tol1 * q), abs(e * q))
            if term1 < term2
                e = d
                d = p / q
            else
                d = xm
                e = d
            end
        else
            d = xm
            e = d
        end

        a = b
        fa = fb

        if (abs(d) > tol1)
            b += d
        else
            b += sign(xm) * abs(tol1)
        end

        fb = f(b)
    end
    return b, fb, iter
end

end # module