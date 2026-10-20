"""
    words_without_k_run(k,n;m=2)

OGF of words without k consecutive occurrences of a designated letter 
for an alphabet of cardinality m (defaults to binary words: m=2).

``\\frac{1 - z^k}{ 1 - mz + (m-1)z^{k+1} }``
For instance, if n=4 and k=3, these words are not counted: aaab, baaa, aaaa.  
"""
function words_without_k_run(k,n;m=2) 
    word_run_ogf(m,k,z) = (1-z^k)/(1 - m*z + (m-1)*z^(k+1))
    taylor_ser = TaylorSeries.taylor_expand(z -> word_run_ogf(m,k,z),0;order=n+1)
    TaylorSeries.getcoeff(taylor_ser,n)
end

"""
    longest_run_binary_asym(n)

Asymptotics for the average of the longest run of ``a``s in a random binary string of length n. 

 ``\\log_2 n`` .
For instance, a random binary word with 10 letters (e.g. baaabababb) will present, 
on average, ``3.32...`` consecutive ``a``s
"""
longest_run_binary_asym(n) = log2(n)


"""
    bin_words_with_k_occurences(k,n)

The set L of binary words that contain exactly k occurrences of the letter b

``L = SEQ(a){(b SEQ(a))}^k  \\implies  L(z) = \\frac{z^k}{{(1-z)}^{k+1}}  ``
For instance, among binary words with 10 letters, there are 210 words  with 4 ``b``s.
(``[z^{10}]L(z) = 210``)
"""
function bin_words_with_k_occurences(k,n)
    word_ogf(z,k) = (z^k)/((1 - z)^(k+1))
    taylor_ser = TaylorSeries.taylor_expand(z->word_ogf(z,k),0;order=n+1)
    TaylorSeries.getcoeff(taylor_ser,n)
end



"""
    bin_words_with_k_occurences_constr(k,n,d)

The set L of binary words that contain exactly k occurrences of the letter b, constrained by the maximum distance d among occurrences

``L^{[d]} = SEQ(a){(b SEQ_{<d}(a))}^{k-1}(b SEQ(a))``.  

For instance, among binary words with 10 letters, there are 45 words with 4 ``b``s in which the maximum distance between ``b``s is 2 (e.g.aaabababab)
(``[z^{10}]L^{[2]}(z) = 45``).
"""
function bin_words_with_k_occurences_constr(k,n,d)
    # Binomial convolution. See page 52 in Flaj  & Sedg.
    accum = Int128[]
    for j in 0:Int(round(n/d))
        if(j>=k-1  || j*d >=n)
            break
        end
        #print("\nj is ",j)
        acc = (-1)^j * binomial(k-1,j) * binomial(n-d*j,k)
        #print("\naccum is ",acc)
        push!(accum,acc)
    end
    #print("\nSum is ")
    sum(accum)
end


"""
    bin_words_runs_coeff(r;n_tot=200)

Taylor series coefficient from generating function for binary words (when p=q=prob(a)=prob(b)) that never have more than r consecutive identical letters. 

The number of binary words that never have more than r consecutive identical letters is found to be (set α = β = r).
n_tot defaults to 200, according to the example in Flajolet & Sedgewick pag. 52
"""
function bin_words_runs_coeff(r;n_tot=200)
    w_rr(r,z) = (1-z^(r+1))/(1-2z+z^(r+1)) # OGF
    #w_rr(r,z) = sum(z^x for x in 0:r)/(1 - sum(z^x for x in 1:r)) # Alternate form

    coefs = TaylorSeries.taylor_expand(z -> w_rr(r,z),0;order=n_tot+1)
    TaylorSeries.getcoeff(coefs,n_tot)
    
end


"""
    bin_words_runs_prob(k,n)

Returns probablity associatied with k-lenght double runs (or either ``a``s or ``b``s) in a sequence of size n (when p=q=prob(a)=prob(b)).   

``W ∼= SEQ(b) SEQ(a SEQ(a) b SEQ(b)) SEQ(a).``
For instance, if n=5 and k=2, the probability of obtaining strings such as bbaab and aabba.  
"""
function bin_words_runs_prob(k,n)
    #a =  FastRational{Int128}(1/(2^n))
    #a = Rational{BigInt}(1/(2^n))
    a = 1//(BigInt(2)^n)
    #a = 1/(2^n)
    result_sym = a*(bin_words_runs_coeff(k;n_tot=n) - bin_words_runs_coeff(k-1;n_tot=n))
    return(Float64(result_sym))
end



"""
    weighted_bin_runs_prob(p,q,l,n)

Weighted model for consecutive runs in binary words.  
    
Probablity of the absence of l-runs among a sequence of n random trials with probabilities p and q.   
``[z^{n}] \\frac{1 - p^l z^l}{1 - z + q p^l z^{l+1}}``.
Use diff over output values to obtain a probability distribution. For n=15, p=0.4 and q=0.6:
`raw_probs = map(x->weighted_bin_runs_prob(0.4,0.6,x,15),0:1:15);plot(diff(raw_probs))`  
"""
function weighted_bin_runs_prob(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    wei_mgf(p,q,l,z) = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1)))
    tay_exp = TaylorSeries.taylor_expand(z -> wei_mgf(p,q,l,z),0;order=n+1)
    TaylorSeries.getcoeff(tay_exp,n)

end



"""
    weighted_bin_runs_pval(p,q,l,n)

p-value obtained from a one-tailed based on the exact distribution using the weighted model for consecutive runs `weighted_bin_runs_coeff`.  
"""
function weighted_bin_runs_pval(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    weighted_coeffs = map(x->weighted_bin_runs_prob(p,q,x,n),(l-1):1:n)
    probs = diff(weighted_coeffs)
    Float64(sum(probs))

end


"""
    weighted_bin_runs_bgf(p, q, l, z, u; mode="count_ge")

Weighted Bivariate Generating Function (BGF) ``W(z, u)`` for consecutive runs in binary words of Bernoulli trials
with probabilities ``p`` and ``q``, where ``z`` marks sequence length ``n`` and ``u`` marks a run parameter
(Flajolet & Sedgewick, Sections I.4 & III.3):

- `mode="count_ge"` (default): ``u`` marks the **number of maximal runs of length ``\\ge l``**:
  ``W(z, u) = \\frac{1 + (u - 1) p^l z^l}{1 - (p + q)z - (u - 1) q p^l z^{l+1}}``
- `mode="length_ge"`: ``u`` marks the **total number of elements belonging to maximal runs of length ``\\ge l``**
  (the numerator of RQA Determinism / Laminarity):
  ``W(z, u) = \\frac{1 + \\beta(z, u)}{1 - qz(1 + \\beta(z, u))}``, where ``\\beta(z, u) = \\frac{pz - p^l z^l}{1 - pz} + \\frac{(upz)^l}{1 - upz}``
- `mode="count_eq"`: ``u`` marks the **number of maximal runs of exact length ``l``**:
  ``W(z, u) = \\frac{1 + (u - 1) p^l z^l (1 - pz)}{1 - (p + q)z - (u - 1) q p^l z^{l+1} (1 - pz)}``
- `mode="double_count_ge"`: ``u`` marks the **number of maximal runs of either outcome (``p`` or ``q``) of length ``\\ge l``**.
"""
function weighted_bin_runs_bgf(p, q, l, z, u; mode="count_ge")
    if mode == "count_ge"
        return (1 + (u - 1) * (p^l) * (z^l)) / (1 - (p + q)*z - (u - 1) * q * (p^l) * (z^(l+1)))
    elseif mode == "length_ge"
        beta_zu = (p*z - (p^l)*(z^l)) / (1 - p*z) + ((u*p*z)^l) / (1 - u*p*z)
        return (1 + beta_zu) / (1 - q*z*(1 + beta_zu))
    elseif mode == "count_eq"
        return (1 + (u - 1) * (p^l) * (z^l) * (1 - p*z)) / (1 - (p + q)*z - (u - 1) * q * (p^l) * (z^(l+1)) * (1 - p*z))
    elseif mode == "double_count_ge"
        alpha_zu = (q*z + (u - 1)*(q^l)*(z^l)) / (1 - q*z)
        beta_zu = (p*z + (u - 1)*(p^l)*(z^l)) / (1 - p*z)
        return (1 + alpha_zu) * (1 + beta_zu) / (1 - alpha_zu * beta_zu)
    else
        error("Unsupported mode '$mode'. Must be one of: \"count_ge\", \"length_ge\", \"count_eq\", \"double_count_ge\".")
    end
end


"""
    weighted_bin_runs_moments(p, q, l, n; mode="count_ge")

Exact analytical moments `(mean, var, std, count)` for run statistics in a weighted binary sequence of length `n`
with probabilities `p` and `q` and run threshold `l`, obtained by differentiating `weighted_bin_runs_bgf` at ``u = 1``.

See [`weighted_bin_runs_bgf`](@ref) for available `mode` options (`"count_ge"`, `"length_ge"`, `"count_eq"`, `"double_count_ge"`).
"""
function weighted_bin_runs_moments(p, q, l, n; mode="count_ge")
    if abs(p + q - 1) > 0.01 || (l > n) || (l < 1)
        println("p + q must be equal to 1 and 1 <= l <= n")
        return (mean = NaN, var = NaN, std = NaN, count = NaN)
    end
    return bgf_moments((z, u) -> weighted_bin_runs_bgf(p, q, l, z, u; mode=mode), n)
end
