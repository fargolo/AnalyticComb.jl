"""
    words_without_k_run(k,n;m=2)

OGF of words without k consecutive occurrences of a designated letter 
for an alphabet of cardinality m (defaults to binary words: m=2).

``\\frac{1 - z^k}{ 1 - mz + (m-1)z^{k+1} }``
For instance, if n=4 and k=3, these words are not counted: aaab, baaa, aaaa.  
"""
function words_without_k_run(k,n;m=2)
    z = SymPy.symbols("z") 
    word_run_ogf = (1-z^k)/(1 - m*z + (m-1)*z^(k+1))
    coefs = collect(series(word_run_ogf,z,0,n+1),z)
    coefs.coeff(z,n)
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
    z = SymPy.symbols("z") 
    word_ogf = (z^k)/((1 - z)^(k+1))
    coefs = collect(series(word_ogf,z,0,n+1),z)
    coefs.coeff(z,n)
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
    W_coeff(r;n_tot=200)

Taylor series coefficient from generating function for binary words that never have more than r consecutive identical letters. 

The number of binary words that never have more than r consecutive identical letters is found to be (set α = β = r).
n_tot defaults to 200, according to the example in Flajolet & Sedgewick pag. 52
"""
function W_coeff(r;n_tot=200)

    z = SymPy.symbols("z") 
    w_rr(r,z) = (1-z^(r+1))/(1-2z+z^(r+1)) # OGF
    #w_rr(r,z) = sum(z^x for x in 0:r)/(1 - sum(z^x for x in 1:r)) # Alternate form

    coefs = collect(series(w_rr(r,z),z,0,n_tot+1),z)
    coefs.coeff(z,n_tot)
    
end


"""
    p_binary_words_doub_runl(k,n)

Returns probablity associatied with k-lenght double runs (or either ``a``s or ``b``s) in a sequence of size n.   

``W ∼= SEQ(b) SEQ(a SEQ(a) b SEQ(b)) SEQ(a).``
For instance, if n=5 and k=2, the probability of obtaining strings such as bbaab and aabba.  
"""
function p_binary_words_doub_runl(k,n)
    #a =  FastRational{Int128}(1/(2^n))
    #a = Rational{BigInt}(1/(2^n))
    a = 1//(BigInt(2)^n)
    #a = 1/(2^n)
    result_sym = a*(W_coeff(k;n_tot=n) - W_coeff(k-1;n_tot=n))
    return(Float64(result_sym))
end




"""
    weighted_bin_runs_coeff(p,q,l,n)

Weighted model for consecutive runs.  
    
Probablity of the absence of l-runs among a sequence of n random trials with probabilities p and q.   
``[z^{n}] \\frac{1 - p^l z^l}{1 - z + q p^l z^{l+1}}``.
Use diff over output values to obtain a probability distribution. For n=15, p=0.4 and q=0.6:
`raw_probs = map(x->weighted_bin_runs_coeff(0.4,0.6,x,15),0:1:15);plot(diff(raw_probs))`  
"""
function weighted_bin_runs_coeff(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    z = SymPy.symbols("z") 
    wei_mgf = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1))) # OGF
    coefs = collect(series(wei_mgf,z,0,n+1),z)
    coefs.coeff(z,n)

end


"""
    p_val_weighted(p,q,l,n)

p-value obtained from a one-tailed based on the exact distribution using the weighted model for consecutive runs `weighted_bin_runs_coeff`.  
"""
function p_val_weighted(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    weighted_coeffs = map(x->weighted_bin_runs_coeff(p,q,x,n),(l-1):1:n)
    probs = diff(weighted_coeffs)
    Float64(sum(probs))

end
