using SymPy
using TaylorSeries
#using FastRationals

@vars z r n

# Binary words that never have more than r consecutive identical letters is found to be (set α = β = r)
# Flajolet pag. 52
w_rr(r,z) = (1-z^(r+1))/(1-2z+z^(r+1)) # OGF
w_rr(r,z) = sum(z^x for x in 0:r)/(1 - sum(z^x for x in 1:r)) # Alternate form

n_words(r;n_tot=200) = getcoeff(taylor_expand(z -> w_rr(z,r),order=n_tot+1),n_tot)
n_words(r;n_tot=200) = getcoeff(taylor_expand(z -> w_rr(z,r),order=n_tot+1),n_tot)


function p_k(k,n)
    #a =  FastRational{Int128}(1/(2^n))
    #a = Rational{BigInt}(1/(2^n))
    a = 1//(BigInt(2)^n)
    #a = 1/(2^n)
    result_sym = a*(n_words(k;n_tot=n) - n_words(k-1;n_tot=n))
    return(Float64(result_sym))
end

# Test
p_k(6,200)
p_k(7,200)

#Flajolet table page 53
for i in 3:12
    print(p_k(i,200))
    print("\n")
end


# p_value com H_0: Probabilidade de série com aquele tamanho apresentar recorrência 
# do tamanho observado ou maior 
# L/Lmax p_k(L,tamanho da coluna)
# V/Vmax (V,tamanho da diagonal)
# TT/NMPRT (TT,tamanho da coluna)


