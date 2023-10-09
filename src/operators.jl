"""
    SEQ(z)

Sequence operator (Pólya quasi-inverse operator).   

Defined as ``A = SEQ(B) \\implies A(z) = \\frac{1}{1 - B(z)}``.
"""
SEQ(z) = 1/(1-z)


"""
    MSET(z)

Multiset operator (Pólya exponential operator).  

Defined as ``A = MSET(B) \\implies A(z) = \\frac{1}{1 - B(z)}``.
"""
function MSET(z,max)
    n = SymPy.symbols("n")
    return(exp(summation(1/n * z^n,(n,1,max))))
end
#MSET(z) = prod(1-z^n)^((-1)*(-B_n))
