# Binary word double runs

Probability for consecutive double runs (either 0s or 1s) of lenght k in binary words of length n, use `p_binary_words_doub_runl(k,n)`

```
julia>using AnalyticComb
julia>p_binary_words_doub_runl(6,200) #e.g. 100000011010... or 01111110101...
0.166...
```

Stirling approximation for n!:
```
julia> stirling_factorial(7)
4980.395831612462
```

# Restricted summands problem

Using SEQ operator to solve Polya partitions with restricted summands (denumerants) problem (see Flajolet & Sedgewick,p. 43) about the number of ways of giving change of 99 cents using pennies (1 cent), nickels (5 cents), dimes (10 cents) and quarters (25 cents).  
```
julia> restricted_sum_part_gf([1,5,10,25]) # examine the generating function SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25)
                 1                  
────────────────────────────────────
        ⎛     5⎞ ⎛     10⎞ ⎛     25⎞
(1 - z)⋅⎝1 - z ⎠⋅⎝1 - z  ⎠⋅⎝1 - z  ⎠

julia>restricted_sum_part(99,[1,5,10,25]) # Counts for 99 as a sum of elements in (1,5,10,25).
213
```


