"""
    SEQ(z)
    SEQ(z, u)

Sequence operator (Pólya quasi-inverse operator).   

Defined as ``A = SEQ(B) \\implies A(z) = \\frac{1}{1 - B(z)}``, or with marker ``u`` for number of components:
``A = SEQ(uB) \\implies A(z, u) = \\frac{1}{1 - u B(z)}`` (Flajolet & Sedgewick, Fig. III.7, p. 170).
"""
SEQ(z) = 1/(1-z)
SEQ(z, u) = 1/(1 - u*z)


"""
    MSET(z,max)
    MSET(z,u,max)

Multiset operator (Pólya exponential operator).  

Defined as ``A = MSET(B) \\implies A(z) = exp(\\sum_{1}^{\\infty} \\frac{1}{k} B(z^k))``, or with marker ``u``
for number of components: ``A(z, u) = exp(\\sum_{1}^{\\infty} \\frac{u^k}{k} B(z^k))``.
Returns a SymPy `:Sym` object.  
"""
function MSET(z,max)
    n = SymPy.symbols("n")
    return(exp(summation(1/n * z^n,(n,1,max))))
end
function MSET(z,u,max)
    n = SymPy.symbols("n")
    return(exp(summation((u^n)/n * z^n,(n,1,max))))
end
#MSET(z) = prod(1-z^n)^((-1)*(-B_n))



"""
    PSET(z,max)
    PSET(z,u,max)

Powerset operator (modified Pólya exponential operator).  

Defined as ``A = PSET(B) \\implies A(z) = exp(\\sum_{1}^{\\infty} \\frac{(-1)^{k-1}}{k} B(z^k))``, or with marker ``u``
for number of components: ``A(z, u) = exp(\\sum_{1}^{\\infty} \\frac{(-1)^{k-1} u^k}{k} B(z^k))``.
Returns a SymPy `:Sym` object.  
"""
function PSET(z,max)
    n = SymPy.symbols("n")
    return(exp(summation((-1)^(n-1)/n * z^n,(n,1,max))))
end
function PSET(z,u,max)
    n = SymPy.symbols("n")
    return(exp(summation((-1)^(n-1)*(u^n)/n * z^n,(n,1,max))))
end


"""
    CYC(z,max)
    CYC(z,u,max)

Cycle operator (Pólya logarithm).  

Defined as ``A = CYC(B) \\implies A(z) = \\sum_{1}^{\\infty} \\frac{\\phi(k)}{k} log \\frac{1}{1-z^k}``, or with marker ``u``
for number of components: ``A(z, u) = \\sum_{1}^{\\infty} \\frac{\\phi(k)}{k} log \\frac{1}{1 - u^k z^k}``.
Returns a SymPy `:Sym` object.  
"""
function CYC(z,max)
    n = SymPy.symbols("n")
    return(summation(SymPy.sympy.totient(n)/n * z^n * log(1/(1-z^n))  ,(n,1,max)))
end
function CYC(z,u,max)
    n = SymPy.symbols("n")
    return(summation(SymPy.sympy.totient(n)/n * log(1/(1 - (u^n)*(z^n))),(n,1,max)))
end


"""
    bgf_moments(bgf::Function, n::Integer)
    bgf_moments(bgf::SymPy.Sym, n::Integer; z=SymPy.symbols("z"), u=SymPy.symbols("u"))

Compute exact analytical moments (`mean`, `var`, `std`, and total weight/count `count`) at size `n`
from a Bivariate Generating Function (BGF) ``A(z, u)`` where ``z`` marks size and ``u`` marks a combinatorial parameter ``\\chi``.

Uses Proposition III.2 of Flajolet & Sedgewick (*Analytic Combinatorics*, p. 158):
- ``\\mathbb{E}_n(\\chi) = \\frac{[z^n] \\partial_u A(z, u)|_{u=1}}{[z^n] A(z, 1)}``
- ``\\mathbb{E}_n(\\chi(\\chi - 1)) = \\frac{[z^n] \\partial_u^2 A(z, u)|_{u=1}}{[z^n] A(z, 1)}``
- ``\\mathbb{V}_n(\\chi) = \\mathbb{E}_n(\\chi(\\chi - 1)) + \\mathbb{E}_n(\\chi) - \\mathbb{E}_n(\\chi)^2``

Returns a `NamedTuple`: `(mean = ..., var = ..., std = ..., count = ...)`.
"""
function bgf_moments(bgf::Function, n::Integer)
    du = TaylorSeries.Taylor1(Float64, 2)
    u = 1.0 + du
    z = TaylorSeries.Taylor1([zero(du), one(du)], Int(n))
    ser = bgf(z, u)
    cn = TaylorSeries.getcoeff(ser, Int(n))
    a_n = Float64(TaylorSeries.getcoeff(cn, 0))
    f1 = Float64(TaylorSeries.getcoeff(cn, 1)) / a_n
    f2 = 2.0 * Float64(TaylorSeries.getcoeff(cn, 2)) / a_n
    mean_val = f1
    var_val = max(0.0, f2 + f1 - f1^2)
    return (mean = mean_val, var = var_val, std = sqrt(var_val), count = a_n)
end

function bgf_moments(bgf::SymPy.Sym, n::Integer; z=SymPy.symbols("z"), u=SymPy.symbols("u"))
    A0_sym = SymPy.subs(bgf, u => 1)
    A1_sym = SymPy.subs(SymPy.diff(bgf, u), u => 1)
    A2_sym = SymPy.subs(SymPy.diff(bgf, u, 2), u => 1)

    f0_fn = SymPy.lambdify(A0_sym, [z])
    f1_fn = SymPy.lambdify(A1_sym, [z])
    f2_fn = SymPy.lambdify(A2_sym, [z])

    z_tay = TaylorSeries.Taylor1(Float64, Int(n))
    a_n = Float64(TaylorSeries.getcoeff(Base.invokelatest(f0_fn, z_tay), Int(n)))
    c1 = Float64(TaylorSeries.getcoeff(Base.invokelatest(f1_fn, z_tay), Int(n)))
    c2 = Float64(TaylorSeries.getcoeff(Base.invokelatest(f2_fn, z_tay), Int(n)))

    f1 = c1 / a_n
    f2 = c2 / a_n
    mean_val = f1
    var_val = max(0.0, f2 + f1 - f1^2)
    return (mean = mean_val, var = var_val, std = sqrt(var_val), count = a_n)
end


# ---------------------------------------------------------------------------
# Phase 1.4: Dominant Singularity Localization & Pole Proximity Metrics
# ---------------------------------------------------------------------------

"""
    dominant_pole(denom_fn::Function; bracket=(0.0, 10.0), n_grid::Int=1000, tol::Real=1e-13) -> Float64

Locate the smallest positive real root ``z_0 > 0`` of a meromorphic/rational generating
function's denominator `denom_fn(z) == 0` via adaptive grid scanning followed by bisection
(Flajolet & Sedgewick, Theorem IV.9 & Section IV.6).

Returns `NaN` if no sign change is detected within `bracket`.
"""
function dominant_pole(denom_fn::Function; bracket=(0.0, 10.0), n_grid::Int=1000, tol::Real=1e-13)
    lo = max(1e-9, Float64(bracket[1]))
    hi = Float64(bracket[2])
    hi > lo || throw(ArgumentError("Upper bracket must exceed lower bracket"))

    grid = collect(range(lo, hi; length=n_grid))
    if lo < 1.0 < hi
        push!(grid, 1.0)
        sort!(grid)
    end

    a = grid[1]
    fa = Float64(denom_fn(a))
    abs(fa) <= tol && return a

    found = false
    b = a
    fb = fa
    for k in 2:length(grid)
        b = grid[k]
        fb = Float64(denom_fn(b))
        if abs(fb) <= tol
            return b
        end
        if fa * fb < 0.0
            found = true
            break
        end
        a = b
        fa = fb
    end

    # Expand upper bound exponentially if denominator is monotonically decaying without crossing yet
    if !found
        step_hi = hi
        prev_x = b
        prev_f = fb
        for _ in 1:20
            step_hi *= 2.0
            f_next = Float64(denom_fn(step_hi))
            if prev_f * f_next <= 0.0
                a, b = prev_x, step_hi
                fa, fb = prev_f, f_next
                found = true
                break
            end
            prev_x, prev_f = step_hi, f_next
        end
    end

    found || return NaN

    for _ in 1:100
        mid = 0.5 * (a + b)
        fmid = Float64(denom_fn(mid))
        if abs(fmid) <= tol || (b - a) <= tol * max(1.0, abs(mid))
            return mid
        end
        if fa * fmid <= 0.0
            b = mid
            fb = fmid
        else
            a = mid
            fa = fmid
        end
    end
    return 0.5 * (a + b)
end


"""
    runs_dominant_pole(p::Real, q::Real, l::Integer)

Compute the dominant real singularity ``z_0`` of the run-avoidance generating function
(irredundant denominator ``1 - q z \\sum_{j=0}^{l-1} (p z)^j = 0``), where ``p + q = 1``
and ``l \\ge 1``.

For run-avoidance, ``P(\\text{no } l\\text{-run in } n) \\sim c \\cdot z_0^{-n}``: as ``p`` increases,
``z_0 > 1`` grows and the avoidance survival rate ``1/z_0`` drops. Conversely, the coherent
persistence pole ``z_{\\text{persist}} = 1/p`` of ``(1 - pz)^{-1}`` satisfies
``z_{\\text{persist}} - 1 \\to 0^+`` as ``p \\to 1`` (Flajolet & Sedgewick, Section V.2–V.3).

Returns `(z0 = ..., gap = z0 - 1.0, growth_rate = 1 / z0, persistence_pole = 1 / p, persistence_gap = 1 / p - 1.0)`.
"""
function runs_dominant_pole(p::Real, q::Real, l::Integer)
    l >= 1 || throw(ArgumentError("Run threshold l must be >= 1, got $l"))
    pf = Float64(p)
    qf = Float64(q)
    abs(pf + qf - 1.0) <= 1e-4 || throw(ArgumentError("p + q must equal 1"))

    pers_pole = pf > 0 ? 1.0 / pf : Inf
    if pf >= 1.0 - 1e-15
        return (z0 = Inf, gap = Inf, growth_rate = 0.0, persistence_pole = 1.0, persistence_gap = 0.0)
    elseif qf >= 1.0 - 1e-15
        return (z0 = 1.0, gap = 0.0, growth_rate = 1.0, persistence_pole = Inf, persistence_gap = Inf)
    end

    denom_fn = function (z)
        s = 0.0
        pz = 1.0
        for _ in 0:(l - 1)
            s += pz
            pz *= pf * z
        end
        return 1.0 - qf * z * s
    end

    z0 = dominant_pole(denom_fn; bracket=(1.0, 20.0))
    return (
        z0 = z0,
        gap = z0 - 1.0,
        growth_rate = 1.0 / z0,
        persistence_pole = pers_pole,
        persistence_gap = pers_pole - 1.0
    )
end


"""
    markov_dominant_pole(p00::Real, p01::Real, p10::Real, p11::Real, l::Integer)

Compute the dominant real singularity ``z_0`` of the 1st-order Markov run-avoidance
generating function (cleared denominator ``(1 - p_{00} z) - p_{01} p_{10} z \\sum_{k=1}^{l-1} p_{11}^{k-1} z^k = 0``),
together with the coherent persistence pole ``z_{\\text{persist}} = 1 / p_{11}`` of ``(1 - p_{11} z)^{-1}``.

As coherent persistence ``p_{11} \\to 1``, the run-avoidance decay rate ``1/z_0`` drops and
``z_{\\text{persist}} - 1 \\to 0^+``.

Returns `(z0 = ..., gap = z0 - 1.0, growth_rate = 1 / z0, persistence_pole = 1 / p11, persistence_gap = 1 / p11 - 1.0)`.
"""
function markov_dominant_pole(p00::Real, p01::Real, p10::Real, p11::Real, l::Integer)
    l >= 1 || throw(ArgumentError("Run threshold l must be >= 1, got $l"))
    p00f = Float64(p00); p01f = Float64(p01)
    p10f = Float64(p10); p11f = Float64(p11)
    pers_pole = p11f > 0 ? 1.0 / p11f : Inf

    if l == 1
        z0 = p00f > 0 ? 1.0 / p00f : Inf
        return (z0 = z0, gap = z0 - 1.0, growth_rate = p00f, persistence_pole = pers_pole, persistence_gap = pers_pole - 1.0)
    end
    if p01f <= 1e-15 || p10f <= 1e-15
        z0 = p00f > 0 ? 1.0 / p00f : Inf
        return (z0 = z0, gap = z0 - 1.0, growth_rate = p00f, persistence_pole = pers_pole, persistence_gap = pers_pole - 1.0)
    end

    denom_fn = function (z)
        s = 0.0
        p11_pow = 1.0
        z_pow = z
        for _ in 1:(l - 1)
            s += p11_pow * z_pow
            p11_pow *= p11f
            z_pow *= z
        end
        return (1.0 - p00f * z) - p01f * p10f * z * s
    end

    z0 = dominant_pole(denom_fn; bracket=(1.0, 20.0))
    return (
        z0 = z0,
        gap = z0 - 1.0,
        growth_rate = 1.0 / z0,
        persistence_pole = pers_pole,
        persistence_gap = pers_pole - 1.0
    )
end
