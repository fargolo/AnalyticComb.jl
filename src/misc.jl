"""
    stirling_factorial_asym(n)

Stirling approximation for n! (EIS A000142)

``n! \\sim \\sqrt{2 \\pi n} {\\frac{n}{e}}^{n}``
"""
stirling_factorial_asym(n) = (n/exp(1))^n*sqrt(2*pi*n)

"""
    stirling_catalan_asym(n)

Stirling approximation for n_th Catalan number. (EIS A000108)

``C_{n} \\sim \\frac{4^n}{\\sqrt{\\pi n^3}}``
"""
stirling_catalan_asym(n) = 4^n/(sqrt(pi*n^3))

"""
    catalan_factorial(n)

The n_th Catalan number. (EIS A000108)

Calculated using factorials.  
``C_{n} = \\frac{(2n)!}{(n+1)!n!}``
"""
catalan_factorial(n) = factorial(2n) / (factorial(n+1) * factorial(n))


"""
    general_trees(n)

The number of general trees, given by (n-1)_th Catalan number. (EIS A000108)

Calculated using binomial coefficients.  
``G_n = C_{n-1} = \\frac{1}{n} {2n - 2 \\choose n - 1}``
"""
general_trees(n) = 1/n * binomial(2n - 2,n-1)


"""
    triangulations(n)

The number of triangulations, given by the n_th Catalan number. (EIS A000108)

Calculated using binomial coefficients. A maximal decomposition of the convex (n + 2)-gon defined by the points into n triangles.    
``T_n = C_{n-1} = \\frac{1}{n+1} {2n \\choose n}``
"""
triangulations(n) = 1/(n+1) * binomial(2n,n)
