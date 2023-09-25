#p 16
binary_count(n) = 2^n

permutation_count(n) = factorial(n)

# p 20 and 35
catalan_gen(z) = (1−sqrt(1−4z))/2z
series(catalan_gen(z),z,0,10)

#p38
stirling_catalan(n) = 4^n/(sqrt(pi*n^3))