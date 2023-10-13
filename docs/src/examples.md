# Binary word double runs

Probability for consecutive double runs (either 0s or 1s) of lenght k in binary words of length n, use `p_binary_words_doub_runl(k,n)`

```
julia>using AnalyticComb
julia>p_binary_words_doub_runl(6,200) #e.g. 100000011010... or 01111110101...
0.166...
```

# Restricted summands problem

Solve Polya partitions with restricted summands (denumerants) problem (Flajolet & Sedgewick, p.43). 
For instance, what is the number of ways of giving change of 99 cents using pennies (1 cent), nickels (5 cents), dimes (10 cents) and quarters (25 cents). 

That is, the number of ways to obtain 99 by summing 1s,5s,10s and 25s.   
``99 = k_1 1 + k_2 5 + k_3 10 + k_4 25 ``.  

The generating function is:``P(z) = SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25)`` and the solution is the
cofficient of ``z^99`` in the expansion: ``[z^99] T(P(z))``.  

Function `restricted_sum_part_gf(r)` returns the generating function for elements in `r` and `restricted_sum_part(k,r)` returns the coefficient in ``z^k`` for the generating function with elements in `r`.  

```
julia> restricted_sum_part_gf([1,5,10,25]) # examine the generating function SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25)
                 1                  
────────────────────────────────────
        ⎛     5⎞ ⎛     10⎞ ⎛     25⎞
(1 - z)⋅⎝1 - z ⎠⋅⎝1 - z  ⎠⋅⎝1 - z  ⎠

julia>restricted_sum_part(99,[1,5,10,25]) # Counts for 99 as a sum of elements in (1,5,10,25).
213
```

