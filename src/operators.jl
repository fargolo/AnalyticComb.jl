"""
    SEQ(z)

Sequence operator (Pólya quasi-inverse operator). 

Defined as ``A = SEQ(B) \\implies A(z) = \\frac{1}{1 - B(z)}``.
"""
function SEQ(z)
    z = SymPy.symbols("z")
    return(1/(1-z))
end


"""
    MSET(z)

Multiset operator (Pólya exponential operator). 

Defined as ``A = MSET(B) \\implies A(z) = \\frac{1}{1 - B(z)}``.
"""
function MSET(z,max)
    z = SymPy.symbols("z")
    n = SymPy.symbols("n")
    return(exp(summation(1/n * z^n,(n,1,max))))
end
#MSET(z) = prod(1-z^n)^((-1)*(-B_n))
