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

Returns probablity associatied with k-lenght double runs (or either 0s or 1s) in a sequence of size n. 

Specification is W ∼= SEQ(b) SEQ(a SEQ(a) b SEQ(b)) SEQ(a).
Refer to the example in Flajolet & Sedgewick pag. 52.  
"""
function p_binary_words_doub_runl(k,n)
    #a =  FastRational{Int128}(1/(2^n))
    #a = Rational{BigInt}(1/(2^n))
    a = 1//(BigInt(2)^n)
    #a = 1/(2^n)
    result_sym = a*(W_coeff(k;n_tot=n) - W_coeff(k-1;n_tot=n))
    return(Float64(result_sym))
end

