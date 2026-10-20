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
