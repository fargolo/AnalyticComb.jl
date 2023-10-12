"""
    SEQ(z)

Sequence operator (Pólya quasi-inverse operator).   

Defined as ``A = SEQ(B) \\implies A(z) = \\frac{1}{1 - B(z)}``.
"""
SEQ(z) = 1/(1-z)


"""
    MSET(z,max)

Multiset operator (Pólya exponential operator).  

Defined as ``A = MSET(B) \\implies A(z) = exp(\\sum{1}_{\\infty} \\frac{1}{k} B(z^k))``.
"""
function MSET(z,max)
    n = SymPy.symbols("n")
    return(exp(summation(1/n * z^n,(n,1,max))))
end
#MSET(z) = prod(1-z^n)^((-1)*(-B_n))



"""
    PSET(z,max)

Powerset operator (modified Pólya exponential operator).  

Defined as ``A = PSET(B) \\implies A(z) = exp(\\sum{1}_{\\infty} \\frac{(-1)^{k-1}}{k} B(z^k))``.
"""
function PSET(z,max)
    n = SymPy.symbols("n")
    return(exp(summation((-1)^(n-1)/n * z^n,(n,1,max))))
end


"""
    CYC(z,max)

Cycle operator (Pólya logarithm).  

Defined as ``A = CYC(B) \\implies A(z) = \\sum{1}_{\\infty} \\frac{\\phi(k)}{k} log \\frac{1}{1-z^k}``.
"""
function CYC(z,max)
    n = SymPy.symbols("n")
    return(summation(SymPy.sympy.totient(n)/n * z^n * log(1/(1-z^n))  ,(n,1,max)))
end
