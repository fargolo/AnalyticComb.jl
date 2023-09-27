# AnalyticComb.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fargolo.github.io/AnalyticComb.jl/dev/)
[![Build Status](https://github.com/fargolo/AnalyticComb.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fargolo/AnalyticComb.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/fargolo/AnalyticComb.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fargolo/AnalyticComb.jl)


# Introduction  
This package implements solutions for combinatorial problems using analytic combinatorics.
Check the text book by Flajojelt & Sedgewick ( https://algo.inria.fr/flajolet/Publications/book.pdf ) and Coursera's full course by Robert Sedgewick ( https://www.coursera.org/learn/analytic-combinatorics ).  

In, 1751, Euler was studying the number of ways in which a given convex polygon could be decomposed into triangles by diagonal lines. (Flajolet & Sedgewick, p.20)

He realized that the progression of numbers in the solution (1, 2, 5, 14, 42, 132,...) was directly related to the coefficients of the series expansion of the polynomial fraction (1−2a−√(1−4a)) / (2aa), that is:
1+2a +5a^2 + 14a^3 + 42a^4 + 132a^5 + ...

Given any constructable combinatorial structure, one can use a set of operators to find a generating function and then approach the problem analytically.

# Install

Python package sympy is required. There's a Julia native implementation in /dev, but currently not working (see  https://discourse.julialang.org/t/trying-julia-for-analytic-combinatorics/103682 ).  

```
$python -m pip install --upgrade pip
$pip install sympy
```

Then, from Julia:  
```
pkg>add AnalyticComb
```

# Quick start  

SymPy.jl functionalities are reexported.  

Probability for consecutive runs of lenght k in binary words of length n, use `p_binary_words_runl(k,n)`

```
julia>using AnalyticComb
julia>p_binary_words_runl(6,200)
0.166...
```

Stirling approximation for n!:
```
julia> stirling_factorial(7)
4980.395831612462
```

Using SEQ operator to solve Polya partitions with restricted summands (denumerants) problem (see Flajolet & Sedgewick,p. 43) about the number of ways of giving change of 99 cents using pennies (1 cent), nickels (5 cents), dimes (10 cents) and quarters (25 cents).  
```
julia> restricted_sum_part_gf([1,5,10,25]) # generating function SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25)
                 1                  
────────────────────────────────────
        ⎛     5⎞ ⎛     10⎞ ⎛     25⎞
(1 - z)⋅⎝1 - z ⎠⋅⎝1 - z  ⎠⋅⎝1 - z  ⎠

julia>restricted_sum_part(99,[1,5,10,25]) 
213
```
