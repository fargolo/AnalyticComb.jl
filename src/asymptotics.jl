"""
    stirling_factorial(n)

Stirling approximation for n! as (n/exp(1))^n*sqrt(2*pi*n).
"""
stirling_factorial(n) = (n/exp(1))^n*sqrt(2*pi*n)

"""
    stirling_catalan(n)

Stirling approximation for n_th Catalan number .
"""
stirling_catalan(n) = 4^n/(sqrt(pi*n^3))

