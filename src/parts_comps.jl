"""
    I_gf(z)

Integers as combinatorial structures 

``I(z)= \\sum_{n \\geq 1} z^n = \\frac{z}{1-z}``.  
Returns a SymPy `:Sym` object.  
"""
function I_gf()
    z = SymPy.symbols("z")
    return(z*SEQ(z))
end


"""
    partitions_gf(z,max)

Generating function for integer partitions.

Defined as ``P(z) = \\prod_{m = 1}^{\\infty} \\frac{1}{1-z^m}``  
Returns a SymPy `:Sym` object.  
Use `SymPy.series` to obtain counts(EIS A000041): `SymPy.series(partitions_gf(z,10),z,0,8)` for n up to 8.  
"""
function partitions_gf(z,max)
    n = SymPy.symbols("n")
    z = SymPy.symbols("z")
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
Returns a SymPy `:Sym` object.      
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
    restricted_sum_comp_gf(z,r) = (1-z)/(1-2z+z^(r+1))
    taylor_ser = TaylorSeries.taylor_expand(z -> restricted_sum_comp_gf(z,r),order=n+1)
    return(TaylorSeries.getcoeff(taylor_ser,n))
end


"""
    restricted_sum_part_gf(r)

Generating function for partition with restricted summand.  
Returns a SymPy `:Sym` object.  
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
    restricted_sum_part_gf(z,r) = prod([SEQ(z^i) for i in r])
    taylor_ser = TaylorSeries.taylor_expand(z -> restricted_sum_part_gf(z,r),order=n+1)
    return(TaylorSeries.getcoeff(taylor_ser,n))
end

"""
    partitions_max_r(n,r)

Number of partitions of size n whose summands lie in the set {1, 2, ... , r} 

n must be an integer and r is the maximum value in the set of summands.
"""
partitions_max_r(n,r) = n^(r-1)*(1 / (factorial(r)*factorial(r-1)))


"""
    partitions_asym(n,tau)

Asymptotics for partitions of size n whose summands lie in the arbitrary finite set ``\\tau`` (tau). 

tau must be an array of integers.
``{P_n}^{T} \\sim \\frac{1}{\\tau} \\frac{n^{r-1}}{(r-1)!}``
"""
partitions_asym(n,tau) = (1/prod(tau))* ((n^(length(tau)-1))/(length(tau) -1))

"""
    fixed_size_comps(n,k)

Compositions of size n made of k summands,

``{C_n}^{k} = {n-1 \\choose k -1}``
"""
fixed_size_comps(n,k) = binomial(n-1,k-1)


"""
    fixed_size_comps_asym(n,k)

Asymptotics for compositions of size n made of k summands,

``{C_n}^{k} \\sim \\frac{n^{k-1}}{(k-1)!}``
"""
fixed_size_comps_asym(n,k) = n^(k-1)/factorial(k-1)