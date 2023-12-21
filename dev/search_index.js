var documenterSearchIndex = {"docs":
[{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"CurrentModule = AnalyticComb","category":"page"},{"location":"functions/#Index","page":"Function definitions","title":"Index","text":"","category":"section"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"List of functions.  ","category":"page"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"","category":"page"},{"location":"functions/#Operators","page":"Function definitions","title":"Operators","text":"","category":"section"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Basic operators SEQ, MSET, PSET, CYC.  ","category":"page"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Modules = [AnalyticComb]\nPages   = [\"operators.jl\"]","category":"page"},{"location":"functions/#AnalyticComb.CYC-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.CYC","text":"CYC(z,max)\n\nCycle operator (Pólya logarithm).  \n\nDefined as A = CYC(B) implies A(z) = sum_1^infty fracphi(k)k log frac11-z^k. Returns a SymPy :Sym object.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.MSET-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.MSET","text":"MSET(z,max)\n\nMultiset operator (Pólya exponential operator).  \n\nDefined as A = MSET(B) implies A(z) = exp(sum_1^infty frac1k B(z^k)). Returns a SymPy :Sym object.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.PSET-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.PSET","text":"PSET(z,max)\n\nPowerset operator (modified Pólya exponential operator).  \n\nDefined as A = PSET(B) implies A(z) = exp(sum_1^infty frac(-1)^k-1k B(z^k)). Returns a SymPy :Sym object.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.SEQ-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.SEQ","text":"SEQ(z)\n\nSequence operator (Pólya quasi-inverse operator).   \n\nDefined as A = SEQ(B) implies A(z) = frac11 - B(z).\n\n\n\n\n\n","category":"method"},{"location":"functions/#Partitions-and-compositions","page":"Function definitions","title":"Partitions and compositions","text":"","category":"section"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Partitions and compositions of integers.  ","category":"page"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Modules = [AnalyticComb]\nPages   = [\"parts_comps.jl\"]","category":"page"},{"location":"functions/#AnalyticComb.I_gf-Tuple{}","page":"Function definitions","title":"AnalyticComb.I_gf","text":"I_gf(z)\n\nIntegers as combinatorial structures \n\nI(z)= sum_n geq 1 z^n = fracz1-z.   Returns a SymPy :Sym object.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.fixed_size_comps-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.fixed_size_comps","text":"fixed_size_comps(n,k)\n\nCompositions of size n made of k summands,\n\nC_n^k = n-1 choose k -1\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.fixed_size_comps_asym-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.fixed_size_comps_asym","text":"fixed_size_comps_asym(n,k)\n\nAsymptotics for compositions of size n made of k summands,\n\nC_n^k sim fracn^k-1(k-1)\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.partitions_asym-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.partitions_asym","text":"partitions_asym(n,tau)\n\nAsymptotics for partitions of size n whose summands lie in the arbitrary finite set tau (tau). \n\ntau must be an array of integers. P_n^T sim frac1tau fracn^r-1(r-1)\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.partitions_asym-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.partitions_asym","text":"partitions_asym(n)\n\nAsymptotics for partition of integers (EIS A000041) by Hardy and Ramanujan, later improved by Rademache\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.partitions_gf-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.partitions_gf","text":"partitions_gf(z,max)\n\nGenerating function for integer partitions.\n\nDefined as P(z) = prod_m = 1^infty frac11-z^m   Returns a SymPy :Sym object.   Use SymPy.series to obtain counts(EIS A000041): SymPy.series(partitions_gf(z,10),z,0,8) for n up to 8.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.partitions_max_r-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.partitions_max_r","text":"partitions_max_r(n,r)\n\nNumber of partitions of size n whose summands lie in the set {1, 2, ... , r} \n\nn must be an integer and r is the maximum value in the set of summands.\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.primes_composition_asym-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.primes_composition_asym","text":"primes_composition_asym(n)\n\nAsymptotics for composition of n into prime parts (EIS A023360).\n\nB_n sim 030365 * 147622^n\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.restricted_sum_comp-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.restricted_sum_comp","text":"restricted_sum_comp(n,r)\n\nNumber of compositions of n with components in the set {1,2,..,r}. \n\nr = 2 yields Fibonnaci numbers (EIS A000045): F_n = F_n-1 + F_n-2.   r>2 yields generalized Fibonacci numbers.\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.restricted_sum_comp_gf-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.restricted_sum_comp_gf","text":"restricted_sum_comp_gf(r)\n\nGenerating function for compositions with restricted summand. Returns a SymPy :Sym object.      \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.restricted_sum_part-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.restricted_sum_part","text":"restricted_sum_part(n,r)\n\nNumber of partitions with components in r with restricted summand n.\n\nn must be an integer and r must be a set of integers, like in r = [1,5,10,25] , n = 99.\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.restricted_sum_part_gf-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.restricted_sum_part_gf","text":"restricted_sum_part_gf(r)\n\nGenerating function for partition with restricted summand.   Returns a SymPy :Sym object.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#Miscellaneous","page":"Function definitions","title":"Miscellaneous","text":"","category":"section"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Miscellaneous functins.    ","category":"page"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Modules = [AnalyticComb]\nPages   = [\"misc.jl\"]","category":"page"},{"location":"functions/#AnalyticComb.catalan_factorial-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.catalan_factorial","text":"catalan_factorial(n)\n\nThe n_th Catalan number. (EIS A000108)\n\nCalculated using factorials.   C_n = frac(2n)(n+1)n\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.general_trees-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.general_trees","text":"general_trees(n)\n\nThe number of general trees, given by (n-1)_th Catalan number. (EIS A000108)\n\nCalculated using binomial coefficients.   G_n = C_n-1 = frac1n 2n - 2 choose n - 1\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.stirling_catalan_asym-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.stirling_catalan_asym","text":"stirling_catalan_asym(n)\n\nStirling approximation for n_th Catalan number. (EIS A000108)\n\nC_n sim frac4^nsqrtpi n^3\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.stirling_factorial_asym-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.stirling_factorial_asym","text":"stirling_factorial_asym(n)\n\nStirling approximation for n! (EIS A000142)\n\nn sim sqrt2 pi n fracne^n\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.triangulations-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.triangulations","text":"triangulations(n)\n\nThe number of triangulations, given by the n_th Catalan number. (EIS A000108)\n\nCalculated using binomial coefficients. A maximal decomposition of the convex (n + 2)-gon defined by the points into n triangles.     T_n = C_n-1 = frac1n+1 2n choose n\n\n\n\n\n\n","category":"method"},{"location":"functions/#Words-and-regular-languages","page":"Function definitions","title":"Words and regular languages","text":"","category":"section"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Words counts, such as number of words with consecutive runs.    ","category":"page"},{"location":"functions/","page":"Function definitions","title":"Function definitions","text":"Modules = [AnalyticComb]\nPages   = [\"words_languages.jl\"]","category":"page"},{"location":"functions/#AnalyticComb.W_coeff-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.W_coeff","text":"W_coeff(r;n_tot=200)\n\nTaylor series coefficient from generating function for binary words that never have more than r consecutive identical letters. \n\nThe number of binary words that never have more than r consecutive identical letters is found to be (set α = β = r). n_tot defaults to 200, according to the example in Flajolet & Sedgewick pag. 52\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.bin_words_with_k_occurences-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.bin_words_with_k_occurences","text":"bin_words_with_k_occurences(k,n)\n\nThe set L of binary words that contain exactly k occurrences of the letter b\n\nL = SEQ(a)(b SEQ(a))^k  implies  L(z) = fracz^k(1-z)^k+1 For instance, among binary words with 10 letters, there are 210 words  with 4 bs. (z^10L(z) = 210)\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.bin_words_with_k_occurences_constr-Tuple{Any, Any, Any}","page":"Function definitions","title":"AnalyticComb.bin_words_with_k_occurences_constr","text":"bin_words_with_k_occurences_constr(k,n,d)\n\nThe set L of binary words that contain exactly k occurrences of the letter b, constrained by the maximum distance d among occurrences\n\nL^d = SEQ(a)(b SEQ_d(a))^k-1(b SEQ(a)).  \n\nFor instance, among binary words with 10 letters, there are 45 words with 4 bs in which the maximum distance between bs is 2 (e.g.aaabababab) (z^10L^2(z) = 45).\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.longest_run_binary_asym-Tuple{Any}","page":"Function definitions","title":"AnalyticComb.longest_run_binary_asym","text":"longest_run_binary_asym(n)\n\nAsymptotics for the average of the longest run of as in a random binary string of length n. \n\nlog_2 n . For instance, a random binary word with 10 letters (e.g. baaabababb) will present,  on average, 332 consecutive as\n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.p_binary_words_doub_runl-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.p_binary_words_doub_runl","text":"p_binary_words_doub_runl(k,n)\n\nReturns probablity associatied with k-lenght double runs (or either as or bs) in a sequence of size n.   \n\nW = SEQ(b) SEQ(a SEQ(a) b SEQ(b)) SEQ(a) For instance, if n=5 and k=2, the probability of obtaining strings such as bbaab and aabba.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.p_val_weighted-NTuple{4, Any}","page":"Function definitions","title":"AnalyticComb.p_val_weighted","text":"p_val_weighted(p,q,l,n)\n\np-value obtained from a one-tailed based on the exact distribution using the weighted model for consecutive runs weighted_bin_runs_coeff.  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.weighted_bin_runs_coeff_nat-NTuple{4, Any}","page":"Function definitions","title":"AnalyticComb.weighted_bin_runs_coeff_nat","text":"weighted_bin_runs_coeff(p,q,l,n)\n\nWeighted model for consecutive runs.  \n\nProbablity of the absence of l-runs among a sequence of n random trials with probabilities p and q.    z^n frac1 - p^l z^l1 - z + q p^l z^l+1. Use diff over output values to obtain a probability distribution. For n=15, p=0.4 and q=0.6: raw_probs = map(x->weighted_bin_runs_coeff(0.4,0.6,x,15),0:1:15);plot(diff(raw_probs))  \n\n\n\n\n\n","category":"method"},{"location":"functions/#AnalyticComb.words_without_k_run-Tuple{Any, Any}","page":"Function definitions","title":"AnalyticComb.words_without_k_run","text":"words_without_k_run(k,n;m=2)\n\nOGF of words without k consecutive occurrences of a designated letter  for an alphabet of cardinality m (defaults to binary words: m=2).\n\nfrac1 - z^k 1 - mz + (m-1)z^k+1  For instance, if n=4 and k=3, these words are not counted: aaab, baaa, aaaa.  \n\n\n\n\n\n","category":"method"},{"location":"examples/#Binary-word-double-runs","page":"Examples","title":"Binary word double runs","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Probability for consecutive double runs (either 0s or 1s) of lenght k in binary words of length n, use p_binary_words_doub_runl(k,n)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia>using AnalyticComb\njulia>p_binary_words_doub_runl(6,200) #e.g. 100000011010... or 01111110101...\n0.166...","category":"page"},{"location":"examples/#Restricted-summands-problem","page":"Examples","title":"Restricted summands problem","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Solving Polya partitions with restricted summands (denumerants) problem[1].  What is the number of ways of giving change of 99 cents using pennies (1 cent), nickels (5 cents), dimes (10 cents) and quarters (25 cents)?","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"That is, the number of ways to obtain 99 by summing 1s,5s,10s and 25s.   ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Or number of distinct sets S = k_1k_2_k_3k_4 quad k in mathbbN satisfying the equation :  ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"1 k_1  + 5 k_2  + 10 k_3  + 25 k_4  = 99.  ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The generating function is:P(z) = SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25) and the solution is the cofficient of z^99 in the expansion: z^99 T(P(z)).  ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"T(P(z)) =  + 213 z^99 + ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Function restricted_sum_part_gf(r) returns the generating function for elements in r and restricted_sum_part(k,r) returns the coefficient in z^k for the generating function with elements in r.  ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> restricted_sum_part_gf([1,5,10,25]) # examine the generating function SEQ(z)*SEQ(z^5)*SEQ(z^10)*SEQ(z^25)\n                 1                  \n────────────────────────────────────\n        ⎛     5⎞ ⎛     10⎞ ⎛     25⎞\n(1 - z)⋅⎝1 - z ⎠⋅⎝1 - z  ⎠⋅⎝1 - z  ⎠\n\njulia>restricted_sum_part(99,[1,5,10,25]) # Counts for 99 as a sum of elements in (1,5,10,25).\n213","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"[1]: Flajolet, P., & Sedgewick, R. (2009). Analytic combinatorics. Cambridge University press. Page 43","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Python package sympy is required for some functionalities. SymPy.jl is reexported in AnalyticComb.jl.   ","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"$python -m pip install --upgrade pip\n$pip install sympy","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Then, from Julia:  ","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"julia> # type the right bracket to enter pkg REPL ']'\npkg>add AnalyticComb","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Or","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"julia>using Pkg; Pkg.add(\"AnalyticComb\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = AnalyticComb","category":"page"},{"location":"#AnalyticComb","page":"Introduction","title":"AnalyticComb","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Documentation for AnalyticComb.","category":"page"},{"location":"#References","page":"Introduction","title":"References","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This package implements solutions for combinatorial problems using analytic combinatorics. Check the text book by Flajolet & Sedgewick and Coursera's full course by Robert Sedgewick.  ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Kudos to Ricardo Bittencourt for his introductory texts on the subject and for helping in an initial implementation.  ","category":"page"},{"location":"#Background","page":"Introduction","title":"Background","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"In 1751, Euler was studying the number of ways in which a given convex polygon could be decomposed into triangles by diagonal lines.[1]","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"He realized that the progression of numbers in the solution (1, 2, 5, 14, 42, 132,...) was directly related to the coefficients of the series expansion of the polynomial fraction   frac12asqrt14a2aa, that is:   1+2a +5a^2 + 14a^3 + 42a^4 + 132a^5 + ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Given any constructable combinatorial structure, one can use a set of operators to find a generating function and then approach the problem analytically.  ","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"For newcomers, this an analytic approach to combinatorial problems. Modelling this type of problem often relies on intuitive arguments. The symbolic method describe such situations with a grammar of operators: Sum, Cartesian product, Sequence, Multiset, Powerset and Cycle. Such operators yield an algebraic expression (e.g. P(z)), called the generating function, which is directly related to the problem via complex analysis. We are generally interested in the coefficients of its series expansion. That is, let the series expansion of P(z) be T(P(z)) = sum_n=1^infty a_n x^n. Then, the values of a_n correspond to the counts of objects of size n in this combinatorial class.  ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For instance, the number of binary words (e.g. abababbabab...) of size n is given by W_n = 2^n. Using the sequence operator (SEQ(A) implies A(z) = frac11-z) , we find the generating function:  ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"W = SEQ(Z+Z) implies W(z) = frac11 - 2z. T(W(z)) = 1 + 2z + 4z^2 + 8z^3 + .  ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This approach can be used to solve complex problems in a systematic way.    ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"[1]: Flajolet, P., & Sedgewick, R. (2009). Analytic combinatorics. Cambridge University press. Page 20  ","category":"page"}]
}
