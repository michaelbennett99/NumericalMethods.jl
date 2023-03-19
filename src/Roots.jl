module Roots

using LinearAlgebra, ForwardDiff

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
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot::Real`: The root of `f`.
- `fxroot::Real`: The value of `f` at the root.
- `xlow::Real`: The lower bound of the root interval.
- `xhigh::Real`: The upper bound of the root interval.
- `niter::Integer`: The number of iterations performed.
"""
function bisect(
        f::Function, xlow::Real, xhigh::Real;
        tol::Real=1.0e-6, maxiter::Integer=30, kwargs...
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
        
        if (diff < tol) 
            break
        end
    end
 
    if (diff > tol)
        error(
            "Did not converge: xlow = $xlow, xhigh = $xhigh, diff = $diff"
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

- `f::Function`: The function to find the root of.
- `x_0::Real`: The first guess.
- `x_1::Real`: The second guess.
- `tol::Real=1e-6`: The tolerance for the root.
- `max_iter::Integer=1000`: The maximum number of iterations to perform.
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
        error("Did not converge: diff = $diff")
    end

    xroot = x_1
    fxroot = f(xroot; kwargs...)

    return xroot, fxroot, iter
end

"""
    func_iter(f, x_0; λ::Real=0, tol=1e-6, max_iter=1000, kwargs...)

Use the functional iteration method to find the root of a function f. This
method works for functions from R^n to R^n, where n ≥ 1.

# Arguments

- `f::Function`: The function to find the root of.
- `x_0::Union{AbstractVector{<:Real}, Real}`: The first guess.
- `λ::Real=0`: The relaxation parameter.
- `tol=1e-6`: The tolerance for the root.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
- `iter`: The number of iterations performed.
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
        error("Did not converge: diff = $diff")
    end

    xroot = x_0
    fxroot = f(xroot; kwargs...)

    return xroot, fxroot, iter
end


"""
    newton(f, x_0; tol=1e-6, max_iter=1000, kwargs...)

Use Newton's method to find the root of a function f: R → R.

# Arguments

- `f::Function`: The function to find the root of.
- `x_0::Real`: The first guess.
- `tol=1e-6`: The tolerance for the root.
- `max_iter=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
- `iter`: The number of iterations performed.
"""
function newton(
        func::Function, x_0::Real;
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
        error("Did not converge: diff = $diff")
    end

    xroot = x_0
    fxroot = f(xroot)

    return xroot, fxroot, iter
end

"""
    newton(f, x_0; tol=1e-6, max_iter=1000, kwargs...)

Use Newton's method to find the root of a multivariate function f: R^n → R.

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

        x_1 = x_0 .- J(x_0) \ f(x_0)

        diff = norm(x_1 - x_0)

        x_0 = x_1

        if iter == max_iter
            break
        end
    end

    if (diff > tol)
        error("Did not converge: diff = $diff")
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
- `x_0::Real`: The first guess.
- `x_1::Real`: The second guess.
- `rtol::Real=1e-6`: The tolerance for the root.
- `ftol::Real=1e-6`: The tolerance for the function value.
- `max_iter::Integer=1000`: The maximum number of iterations to perform.
- `kwargs...`: Additional keyword arguments to pass to `f`.

# Returns

- `xroot`: The root of `f`.
- `fxroot`: The value of `f` at the root.
- `iter`: The number of iterations performed.
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