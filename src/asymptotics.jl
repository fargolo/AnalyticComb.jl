"""
    stirling_factorial(n)

Stirling approximation for n! as (n/exp(1))^n*sqrt(2*pi*n).
"""
stirling_factorial(n) = (n/exp(1))^n*sqrt(2*pi*n)


"""
    partitions(n)

Asymptotics for partition of integers by Hardy and Ramanujan, later improved by Rademache
"""
partitions(n) = (1/(4*n*sqrt(3)))*exp(pi*sqrt(2n/3))


"""
    primes_composition(n)

Asymptotics for composition of primes.
"""
primes_composition(n) = 0.30365 * 1.47622^n