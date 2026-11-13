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


# ---------------------------------------------------------------------------
# Phase 1.1: 1st-Order Markov Binary Run BGFs, Probabilities & Moments
# ---------------------------------------------------------------------------

"""
    estimate_markov_params(seq::AbstractVector{<:Integer}; laplace::Real=0.0)
    estimate_markov_params(seqs::AbstractVector{<:AbstractVector{<:Integer}}; laplace::Real=0.0)

Estimate 1st-order Markov transition probabilities `(p00, p01, p10, p11, p_stat, q_stat)`
from a binary sequence (or collection of sequences) with elements in `{0, 1}`.
Optional `laplace` adds pseudocount smoothing to each transition cell.
"""
function estimate_markov_params(seq::AbstractVector{<:Integer}; laplace::Real=0.0)
    return estimate_markov_params([seq]; laplace=laplace)
end

function estimate_markov_params(seqs::AbstractVector{<:AbstractVector{<:Integer}}; laplace::Real=0.0)
    n00 = Float64(laplace)
    n01 = Float64(laplace)
    n10 = Float64(laplace)
    n11 = Float64(laplace)
    total_ones = 0
    total_len = 0

    for s in seqs
        isempty(s) && continue
        for x in s
            (x == 0 || x == 1) || throw(ArgumentError("Sequence elements must be 0 or 1, got $x"))
        end
        total_ones += count(==(1), s)
        total_len += length(s)
        for i in 1:(length(s) - 1)
            a, b = s[i], s[i + 1]
            if a == 0 && b == 0
                n00 += 1.0
            elseif a == 0 && b == 1
                n01 += 1.0
            elseif a == 1 && b == 0
                n10 += 1.0
            else
                n11 += 1.0
            end
        end
    end

    emp_p = total_len > 0 ? total_ones / total_len : 0.5
    row0 = n00 + n01
    row1 = n10 + n11

    p01 = row0 > 0 ? n01 / row0 : emp_p
    p00 = 1.0 - p01
    p10 = row1 > 0 ? n10 / row1 : (1.0 - emp_p)
    p11 = 1.0 - p10

    denom = p01 + p10
    p_stat = denom > 0 ? p01 / denom : emp_p
    return (p00 = p00, p01 = p01, p10 = p10, p11 = p11, p_stat = p_stat, q_stat = 1.0 - p_stat)
end


"""
    markov_bin_runs_bgf(p00, p01, p10, p11, l, z, u; mode="count_ge", init_probs=nothing)

Exact 1st-order Markov Bivariate Generating Function ``W_{\\text{Markov}}(z, u)`` for consecutive
runs in binary sequences, where ``z`` marks sequence length ``n`` and ``u`` marks the run statistic
(Flajolet & Sedgewick, Sections I.4.2 & III.6.1):

- `mode="count_ge"` (default): number of maximal 1-runs of length ``\\ge l``
- `mode="length_ge"`: total number of elements belonging to 1-runs of length ``\\ge l``
- `mode="count_eq"`: number of maximal 1-runs of exact length ``l``
- `mode="double_count_ge"`: number of maximal runs of either 0 or 1 of length ``\\ge l``

`init_probs` may be `(q0, p1)` giving initial probabilities ``(P(X_1=0), P(X_1=1))``;
defaults to the stationary distribution ``q_{\\text{stat}} = 1 - p_{\\text{stat}}``,
``p_{\\text{stat}} = p_{01} / (p_{01} + p_{10})``.
"""
function markov_bin_runs_bgf(p00, p01, p10, p11, l::Integer, z, u; mode="count_ge", init_probs=nothing)
    l >= 1 || throw(ArgumentError("Run threshold l must be >= 1, got $l"))
    if init_probs === nothing
        denom = p01 + p10
        p1 = denom > 0 ? p01 / denom : 0.5
        q0 = 1.0 - p1
    else
        q0 = init_probs[1]
        p1 = init_probs[2]
    end

    p11_pow = l == 1 ? one(p11) : p11^(l - 1)
    p00_pow = l == 1 ? one(p00) : p00^(l - 1)
    zl = z^l

    if mode == "count_ge"
        Z_block = z / (1 - p00 * z)
        O_block = ((z - p11_pow * zl) + u * p11_pow * zl) / (1 - p11 * z)
    elseif mode == "length_ge"
        Z_block = z / (1 - p00 * z)
        O_short = (z - p11_pow * zl) / (1 - p11 * z)
        O_long = (u^l) * p11_pow * zl / (1 - u * p11 * z)
        O_block = O_short + O_long
    elseif mode == "count_eq"
        Z_block = z / (1 - p00 * z)
        O_block = z / (1 - p11 * z) + (u - 1) * p11_pow * zl
    elseif mode == "double_count_ge"
        Z_block = ((z - p00_pow * zl) + u * p00_pow * zl) / (1 - p00 * z)
        O_block = ((z - p11_pow * zl) + u * p11_pow * zl) / (1 - p11 * z)
    else
        error("Unsupported mode '$mode'. Must be one of: \"count_ge\", \"length_ge\", \"count_eq\", \"double_count_ge\".")
    end

    core = 1 / (1 - p01 * p10 * Z_block * O_block)
    return 1 + (q0 * Z_block + p1 * O_block + (q0 * p01 + p1 * p10) * Z_block * O_block) * core
end


function _validate_markov_probs(p00, p01, p10, p11)
    if abs(p00 + p01 - 1.0) > 1e-4 || abs(p10 + p11 - 1.0) > 1e-4 ||
       p00 < -1e-9 || p01 < -1e-9 || p10 < -1e-9 || p11 < -1e-9
        throw(ArgumentError("Invalid Markov transition probabilities: rows (p00+p01, p10+p11) must be non-negative and sum to 1"))
    end
end


"""
    markov_bin_runs_moments(p00, p01, p10, p11, l, n; mode="count_ge", init_probs=nothing)

Exact analytical moments `(mean, var, std, count)` for run statistics of threshold `l`
in a 1st-order Markov binary sequence of length `n`.
"""
function markov_bin_runs_moments(p00, p01, p10, p11, l::Integer, n::Integer; mode="count_ge", init_probs=nothing)
    _validate_markov_probs(p00, p01, p10, p11)
    l >= 1 || throw(ArgumentError("Run threshold l must be >= 1, got $l"))
    n >= 0 || throw(ArgumentError("Sequence length n must be >= 0, got $n"))
    if n == 0 || (l > n)
        return (mean = 0.0, var = 0.0, std = 0.0, count = 1.0)
    end
    return bgf_moments(
        (z, u) -> markov_bin_runs_bgf(p00, p01, p10, p11, l, z, u; mode=mode, init_probs=init_probs),
        n
    )
end


"""
    markov_bin_runs_prob(p00, p01, p10, p11, l, n; init_probs=nothing) -> Float64

Exact probability of observing **zero** 1-runs of length ``\\ge l`` in a 1st-order Markov
binary sequence of length `n` (i.e., ``P(L_n < l)`` where ``L_n`` is the longest 1-run).
"""
function markov_bin_runs_prob(p00, p01, p10, p11, l::Integer, n::Integer; init_probs=nothing)
    _validate_markov_probs(p00, p01, p10, p11)
    n >= 0 || throw(ArgumentError("Sequence length n must be >= 0, got $n"))
    l <= 0 && return 0.0
    (n == 0 || l > n) && return 1.0
    z_tay = TaylorSeries.Taylor1(Float64, Int(n))
    gf = markov_bin_runs_bgf(p00, p01, p10, p11, l, z_tay, 0.0; mode="count_ge", init_probs=init_probs)
    return clamp(Float64(TaylorSeries.getcoeff(gf, Int(n))), 0.0, 1.0)
end


"""
    markov_bin_runs_pval(p00, p01, p10, p11, l, n; init_probs=nothing) -> Float64

Exact one-sided upper-tail ``p``-value ``P(L_n \\ge l) = 1 - P(L_n < l)`` for observing
at least one consecutive 1-run of length ``\\ge l`` in a 1st-order Markov binary sequence of length `n`.
"""
function markov_bin_runs_pval(p00, p01, p10, p11, l::Integer, n::Integer; init_probs=nothing)
    _validate_markov_probs(p00, p01, p10, p11)
    n >= 0 || throw(ArgumentError("Sequence length n must be >= 0, got $n"))
    l <= 0 && return 1.0
    l > n && return 0.0
    return clamp(1.0 - markov_bin_runs_prob(p00, p01, p10, p11, l, n; init_probs=init_probs), 0.0, 1.0)
end


# ---------------------------------------------------------------------------
# Phase 1.2: Guibas-Odlyzko Autocorrelation Polynomials & Pattern BGFs
# ---------------------------------------------------------------------------

function _letter_prob(probs::AbstractDict, sym)
    return probs[sym]
end

function _letter_prob(probs::AbstractVector, sym::Integer)
    return probs[sym]
end

"""
    autocorrelation_poly(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, z)

Compute the weighted Guibas–Odlyzko autocorrelation (correlation) polynomial ``\\gamma_{\\mathfrak{p}}(z)``
for an arbitrary pattern ``\\mathfrak{p} = p_1 \\dots p_k`` under letter weights `probs`
(Flajolet & Sedgewick, Section I.4.1, p. 60 & Note III.39, p. 213):

``\\gamma_{\\mathfrak{p}}(z) = \\sum_{j : p_1 \\dots p_{k-j} = p_{j+1} \\dots p_k} \\left(\\prod_{t=k-j+1}^k \\pi_{p_t}\\right) z^j``.
"""
function autocorrelation_poly(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, z)
    k = length(pattern)
    k >= 1 || throw(ArgumentError("Pattern must be non-empty"))
    acc = one(z)
    for j in 1:(k - 1)
        if @view(pattern[1:(k - j)]) == @view(pattern[(j + 1):k])
            w = prod(_letter_prob(probs, pattern[t]) for t in (k - j + 1):k)
            acc = acc + w * (z^j)
        end
    end
    return acc
end


"""
    pattern_occurrences_bgf(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, z, u; overlapping::Bool=true)

Bivariate Generating Function ``W_{\\mathfrak{p}}(z, u)`` for an arbitrary pattern ``\\mathfrak{p}``
of length ``k`` in Bernoulli sequences of symbols with probabilities `probs`, where ``z`` marks
sequence length and ``u`` marks the number of occurrences of ``\\mathfrak{p}``
(Flajolet & Sedgewick, Note III.39, p. 213):

- `overlapping=true` (default): counts possibly overlapping occurrences:
  ``W(z, u) = \\frac{(u - 1) \\gamma(z) - u}{(1 - z)((u - 1)\\gamma(z) - u) + (u - 1)\\pi[\\mathfrak{p}] z^k}``
- `overlapping=false`: counts greedy non-overlapping (renewal) occurrences:
  ``W^{\\text{non-ov}}(z, u) = \\frac{\\gamma(z)}{(1 - z)\\gamma(z) + (1 - u)\\pi[\\mathfrak{p}] z^k}``
"""
function pattern_occurrences_bgf(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, z, u; overlapping::Bool=true)
    k = length(pattern)
    k >= 1 || throw(ArgumentError("Pattern must be non-empty"))
    gamma = autocorrelation_poly(pattern, probs, z)
    pi_p = prod(_letter_prob(probs, s) for s in pattern)
    zk = z^k
    if overlapping
        num = (u - 1) * gamma - u
        den = (1 - z) * num + (u - 1) * pi_p * zk
        return num / den
    else
        return gamma / ((1 - z) * gamma + (1 - u) * pi_p * zk)
    end
end


"""
    pattern_occurrences_moments(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, n::Integer; overlapping::Bool=true)

Exact analytical moments `(mean, var, std, count)` for the number of occurrences of arbitrary
`pattern` in a random Bernoulli word of length `n` with letter probabilities `probs`.
"""
function pattern_occurrences_moments(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, n::Integer; overlapping::Bool=true)
    n >= 0 || throw(ArgumentError("Sequence length n must be >= 0, got $n"))
    if n < length(pattern)
        return (mean = 0.0, var = 0.0, std = 0.0, count = 1.0)
    end
    return bgf_moments(
        (z, u) -> pattern_occurrences_bgf(pattern, probs, z, u; overlapping=overlapping),
        n
    )
end


"""
    pattern_waiting_prob(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, n::Integer) -> Float64

Exact probability that `pattern` does **not** occur anywhere in a random Bernoulli sequence
of length `n` (``[z^n] W_{\\mathfrak{p}}(z, 0)``).
"""
function pattern_waiting_prob(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, n::Integer)
    n >= 0 || throw(ArgumentError("Sequence length n must be >= 0, got $n"))
    n < length(pattern) && return 1.0
    z_tay = TaylorSeries.Taylor1(Float64, Int(n))
    gf = pattern_occurrences_bgf(pattern, probs, z_tay, 0.0; overlapping=true)
    return clamp(Float64(TaylorSeries.getcoeff(gf, Int(n))), 0.0, 1.0)
end


"""
    pattern_pval(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, k_obs::Integer, n::Integer; overlapping::Bool=true) -> Float64

Exact combinatorial upper-tail ``p``-value ``P(X_n \\ge k_{\\text{obs}})`` for observing at least
`k_obs` occurrences of `pattern` in a random Bernoulli sequence of length `n`, computed directly
by extracting ``[z^n u^s] W_{\\mathfrak{p}}(z, u)`` for ``s = 0, \\dots, k_{\\text{obs}}-1``.
"""
function pattern_pval(pattern::AbstractVector, probs::Union{AbstractDict,AbstractVector}, k_obs::Integer, n::Integer; overlapping::Bool=true)
    n >= 0 || throw(ArgumentError("Sequence length n must be >= 0, got $n"))
    k_obs <= 0 && return 1.0
    k = length(pattern)
    max_occ = overlapping ? max(0, Int(n) - k + 1) : div(Int(n), k)
    k_obs > max_occ && return 0.0
    if k_obs == 1
        return clamp(1.0 - pattern_waiting_prob(pattern, probs, n), 0.0, 1.0)
    end
    u_var = TaylorSeries.Taylor1(Float64, Int(k_obs - 1))
    z_var = TaylorSeries.Taylor1([zero(u_var), one(u_var)], Int(n))
    ser = pattern_occurrences_bgf(pattern, probs, z_var, u_var; overlapping=overlapping)
    cn = TaylorSeries.getcoeff(ser, Int(n))
    cdf_below = sum(Float64(TaylorSeries.getcoeff(cn, s)) for s in 0:(Int(k_obs) - 1))
    return clamp(1.0 - cdf_below, 0.0, 1.0)
end
