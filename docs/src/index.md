```@meta
CurrentModule = AnalyticComb
```

# AnalyticComb

Documentation for [AnalyticComb](https://github.com/fargolo/AnalyticComb.jl).

# References   

This package implements solutions for combinatorial problems using analytic combinatorics.
Check the text book by Flajojelt & Sedgewick ( https://algo.inria.fr/flajolet/Publications/book.pdf ) and Coursera's full course by Robert Sedgewick ( https://www.coursera.org/learn/analytic-combinatorics ).  

I thank Ricardo Bittencourt ( https://github.com/ricbit/ ) for his introductory texts on the subject and for helping in an initial implementation.  

# Background  

In 1751, Euler was studying the number of ways in which a given convex polygon could be decomposed into triangles by diagonal lines. (Flajolet & Sedgewick, p.20)

He realized that the progression of numbers in the solution (1, 2, 5, 14, 42, 132,...) was directly related to the coefficients of the series expansion of the polynomial fraction 
``\frac{1−2a−\sqrt{1−4a}}{2aa}``, that is:
``1+2a +5a^2 + 14a^3 + 42a^4 + 132a^5 + ...``

Given any constructable combinatorial structure, one can use a set of operators to find a generating function and then approach the problem analytically.

# Introduction  

For newcomers, this an analytic approach to combinatorial problems. Modelling this type of problem often relies on intuitive arguments. Analytic combinatorics describe such situations with a grammar of operators: Sum, Cartesian product, Sequence, Multiset, Powerset and Cycle. Such operators yield an algebraic expression (e.g. ``P(z)``), called the generating function, which is directly related to the problem via complex analysis. We are generally interested in the coefficients of its series expansion. That is, let the series expansion of ``P(z)`` be ``T(P(z)) = \sum_{n=1}^{\infty} a_n x^n``. Then, the values of ``a^n`` correspond to the counts of objects of size ``n`` in this combinatorial class.

For instance, the number of binary words (e.g. abababbabab...) of size n is given by ``W_n = 2^n``. Using the sequence operator (``SEQ(A) \implies A(z) = \frac{1}{1-z}``) , we find the generating function: 

``W = SEQ(Z+Z) \implies W(z) = \frac{1}{1 - 2z}``. ``T(W(z)) = 1 + 2z + 4z^2 + 8z^3 + ...``.  

This approach can be used to solve complex problems in a systematic way.    


