"""
    I_gf(z)

Integers as combinatorial structures 

``I(z)= \\sum_{n \\geq 1} z^n = \\frac{z}{1-z}``
"""
function I_gf(z)
    z = SymPy.symbols("z")
    return(z*SEQ(z))
end

"""
    partitions_gf(z,max)

Generating function for integer partitions.

Defined as ``P(z) = \\prod{m = 1}_{\\Inf} \\frac{1}{1-z^m}``  
Use `series` to obtain counts(EIS A000041): `series(partitions_gf(z,10),z,0,8)` for n up to 8.
"""
function partitions_gf(z,max)
    n = SymPy.symbols("n")
    prod([1/(1-z^n) for n in 1:max])
end
#MSET(I_gf(z)) yields the correct gen. function, as in page 41 (equation 38), 
# but series(MSET(I_gf(z))) returns coefficients in 2^n


"""
    partitions_asym(n)

Asymptotics for partition of integers (EIS A000041) by Hardy and Ramanujan, later improved by Rademache
"""
partitions_asym(n) = (1/(4*n*sqrt(3)))*exp(pi*sqrt(2n/3))


"""
    primes_composition_asym(n)

Asymptotics for composition of n into prime parts (EIS A023360).

``B_{n} \\sim 0.30365 * 1.47622^n``
"""
primes_composition_asym(n) = 0.30365 * 1.47622^n


"""
    restricted_sum_comp_gf(r)

Generating function for compositions with restricted summand.
"""
function restricted_sum_comp_gf(r)
    z = SymPy.symbols("z")
    return((1-z)/(1-2z+z^(r+1)))
end

"""
    restricted_sum_comp(n,r)

Number of compositions of n with components in the set {1,2,..,r}. 

r = 2 yields Fibonnaci numbers (EIS A000045): ``F_{n} = F_{n-1} + F_{n-2}``.  
r>2 yields generalized Fibonacci numbers.
"""
function restricted_sum_comp(n,r)
    z = SymPy.symbols("z")
    coefs = collect(series(restricted_sum_comp_gf(r),z,0,n+1),z)
    return(coefs.coeff(z,n))
end


"""
    restricted_sum_part_gf(r)

Generating function for partition with restricted summand.
"""
function restricted_sum_part_gf(r)
    z = SymPy.symbols("z")
    return(prod([SEQ(z^i) for i in r]))
end

"""
    restricted_sum_part(n,r)

Number of partitions with components in r with restricted summand n.

n must be an integer and r must be a set of integers, like in r = [1,5,10,25] , n = 99.
"""
function restricted_sum_part(n,r)
    z = SymPy.symbols("z")
    coefs = collect(series(restricted_sum_part_gf(r),z,0,n+1),z)
    return(coefs.coeff(z,n))
end