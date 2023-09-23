# p27

sum(a,b) = a + b
cartesian_prod(a,b) = a*b
SEQ(z) = 1/(1-z)
PSET(b(z)) = 
MSET(b(z)) =
CYC(b(z)) = 

#p 34
binary_words_ogf = SEQ(2z)
series(binary_words_ogf)

general_trees_ogf(z) = z*SEQ(general_trees_ogf(z))
general_trees_ogf = g*(1-g) - z 
solve(general_trees_ogf,g) # Solutions are Catalan OGF
series(solve(general_trees_ogf,g)[1],z,0,10) # Catalan numebers

#p36
triangulations_ogf(z) = 1/2z*(1 - sqrt(1 -4z))
series(triangulations_ogf(z)) # Also Catalan numbers

# p37
#Stirling
stirling_aprox(n) = (n/exp(1))^n*sqrt(2*pi*n)

#41
# Compositions
compositions(n) = 2^(n-1)
# Partitions by Hardy and Ramanujan later improved by Rademache
partitions(n) = (1/(4*n*sqrt(3)))*exp(pi*sqrt(2n/3))

#p43
primes_composition(n) = 0.30365 * 1.47622^n

#change
using SymPy
change_counts(z) = 1/((1-z)(1 - z^5)(1 - z^10)(1 - z^25))
series(change_counts(z),z,0,99)
coefs = collect(series(change_counts(z),z,0,99),z)
coefs.coeff(z,99)
coefs.coeff(z^99)
coefs.coeffs()
@vars z
p.coeffs()
