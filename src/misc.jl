"""
    stirling_factorial_asym(n)

Stirling approximation for n! as (n/exp(1))^n*sqrt(2*pi*n).

``n! \\sim \\sqrt{2 \\pi n} {\\frac{n}{e}}^{n}``
"""
stirling_factorial_asym(n) = (n/exp(1))^n*sqrt(2*pi*n)

"""
    stirling_catalan_asym(n)

Stirling approximation for n_th Catalan number .

``frac{4^n}{\\sqrt{\\pi n^3}}``
"""
stirling_catalan_asym(n) = 4^n/(sqrt(pi*n^3))

