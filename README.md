# AnalyticComb.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fargolo.github.io/AnalyticComb.jl/dev/)
[![Build Status](https://github.com/fargolo/AnalyticComb.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fargolo/AnalyticComb.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/fargolo/AnalyticComb.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fargolo/AnalyticComb.jl)


# Introduction  

In 1751, Euler was studying the number of ways in which a given convex polygon could be decomposed into triangles by diagonal lines.[^1]  

He realized that the progression of numbers in the solution $S = 1, 2, 5, 14, 42, 132,...$ was directly related to the coefficients of the series expansion of the polynomial fraction $\frac{1 − 2a − \sqrt{1−4a}}{2aa}$, that is:
$1 + 2a + 5a^2 + 14a^3 + 42a^4 + 132a^5 + ...$

Given any constructable combinatorial structure, one can use a set of operators to find a generating function and then approach the problem analytically.

See the [docs](https://fargolo.github.io/AnalyticComb.jl/dev/).

Check the text book by Flajolet & Sedgewick ( https://algo.inria.fr/flajolet/Publications/book.pdf ) and Coursera's full course by Robert Sedgewick ( https://www.coursera.org/learn/analytic-combinatorics ) for more.  

Kudos to Ricardo Bittencourt ( https://github.com/ricbit/ ) for his introductory texts on the subject and for helping in an initial implementation.

# Quick start  

## Install

Python package sympy is required. 

```
$python -m pip install --upgrade pip
$pip install sympy
```

Then, from Julia:  
```
pkg>add AnalyticComb
```

## Example

This software can be used to solve problems such as Polya's problem of partitions with restricted summands [^2]. What is the number of ways of giving change of 99 cents using pennies (1 cent), nickels (5 cents), dimes (10 cents) and quarters (25 cents)?

```
julia> using AnalyticComb
julia> restricted_sum_part_gf([1,5,10,25]) # examine the generating function from specification SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25)
                 1                  
────────────────────────────────────
        ⎛     5⎞ ⎛     10⎞ ⎛     25⎞
(1 - z)⋅⎝1 - z ⎠⋅⎝1 - z  ⎠⋅⎝1 - z  ⎠

julia>restricted_sum_part(99,[1,5,10,25]) # Counts for 99 as a sum of elements in (1,5,10,25).
213
```

[^1]: Flajolet, P., & Sedgewick, R. (2009). Analytic combinatorics. Cambridge University press. Page 20
[^2]: Ibid. page 43
