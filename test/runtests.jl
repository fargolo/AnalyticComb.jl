using AnalyticComb
using Test

@testset "AnalyticComb.jl" begin

    @testset "Words" begin
       
        # Binary words with distance constraints
        ## 4 occurences in a 10-letter word separated by maximum distance of 2 
        @test bin_words_with_k_occurences_constr(4,10,2) == 45

        @test words_without_k_run(3,10) == 504

        binary_double_results = map(i->bin_words_runs_prob(i,200),3:12)
             
        # Check if any of the results presents error > 0.01
        @test sum(binary_double_results - map(x-> round(x;digits=3),
            [6.5494e-8,0.00070,0.03397,0.16602,0.25741,0.22352,0.14594,0.08292,0.04402,0.02260]) .< 0.01) != 0
        
    end
    
    @testset "PartComps" begin

        
        @test restricted_sum_part(99,[1,5,10,25]) == 213 

        z = SymPy.symbols("z")

        integers = collect(SymPy.series(I_gf(),z,0,10),z) 
        @test integers.coeff(z,3) == 1

        A000041 = [1, 2, 3, 5, 7, 11, 15, 22, 30, 42]
        part_asym_diffs = (map(x -> partitions_asym(x),collect(1:10)) - A000041) ./ A000041 
        @test !(0 in (part_asym_diffs .< 1))

        @test fixed_size_comps(10,3) == 36
        @test fixed_size_comps_asym(10,3) == 50


    end

    @testset "BGF Moments and Marked Runs" begin
        # Example III.5 (Flajolet & Sedgewick, p. 167): Number of summands in integer compositions
        # C(z, u) = (1 - z) / (1 - z*(u + 1)); for n=10: count=2^9=512, mean=(n+1)/2=5.5, var=(n-1)/4=2.25, std=1.5
        comp_mom = bgf_moments((z, u) -> (1 - z) / (1 - z * (u + 1)), 10)
        @test isapprox(comp_mom.count, 512.0; atol=1e-9)
        @test isapprox(comp_mom.mean, 5.5; atol=1e-9)
        @test isapprox(comp_mom.var, 2.25; atol=1e-9)
        @test isapprox(comp_mom.std, 1.5; atol=1e-9)

        # SymPy overload of bgf_moments on marked SEQ(I_gf(), u)
        z, u = SymPy.symbols("z u")
        comp_sym_bgf = SEQ(z / (1 - z), u)
        comp_sym_mom = bgf_moments(comp_sym_bgf, 10; z=z, u=u)
        @test isapprox(comp_sym_mom.mean, 5.5; atol=1e-9)
        @test isapprox(comp_sym_mom.var, 2.25; atol=1e-9)

        # Verify weighted_bin_runs_bgf at u=0 recovers weighted_bin_runs_prob (probability of 0 runs >= l)
        p, q, l, n = 0.4, 0.6, 3, 8
        z_tay = AnalyticComb.TaylorSeries.Taylor1(Float64, n)
        prob_zero_runs = AnalyticComb.TaylorSeries.getcoeff(weighted_bin_runs_bgf(p, q, l, z_tay, 0.0; mode="count_ge"), n)
        @test isapprox(prob_zero_runs, weighted_bin_runs_prob(p, q, l, n); atol=1e-12)

        # Brute-force verification of all 4 modes of weighted_bin_runs_moments against exact enumeration over {0,1}^n
        exact_count_ge = 0.0
        exact_count_ge_sq = 0.0
        exact_len_ge = 0.0
        exact_len_ge_sq = 0.0
        exact_count_eq = 0.0
        exact_count_eq_sq = 0.0
        exact_double_ge = 0.0
        exact_double_ge_sq = 0.0

        for mask in 0:(2^n - 1)
            bits = [(mask >> k) & 1 for k in 0:(n - 1)]
            k1 = count(==(1), bits)
            w = (p^k1) * (q^(n - k1))

            # Extract maximal runs
            runs_1 = Int[]
            runs_all = Int[]
            cur_val = bits[1]
            cur_len = 1
            for idx in 2:n
                if bits[idx] == cur_val
                    cur_len += 1
                else
                    push!(runs_all, cur_len)
                    cur_val == 1 && push!(runs_1, cur_len)
                    cur_val = bits[idx]
                    cur_len = 1
                end
            end
            push!(runs_all, cur_len)
            cur_val == 1 && push!(runs_1, cur_len)

            c_ge = count(r -> r >= l, runs_1)
            l_ge = sum(r for r in runs_1 if r >= l; init=0)
            c_eq = count(r -> r == l, runs_1)
            d_ge = count(r -> r >= l, runs_all)

            exact_count_ge += w * c_ge
            exact_count_ge_sq += w * c_ge^2
            exact_len_ge += w * l_ge
            exact_len_ge_sq += w * l_ge^2
            exact_count_eq += w * c_eq
            exact_count_eq_sq += w * c_eq^2
            exact_double_ge += w * d_ge
            exact_double_ge_sq += w * d_ge^2
        end

        m_cge = weighted_bin_runs_moments(p, q, l, n; mode="count_ge")
        @test isapprox(m_cge.mean, exact_count_ge; atol=1e-10)
        @test isapprox(m_cge.var, exact_count_ge_sq - exact_count_ge^2; atol=1e-10)

        m_lge = weighted_bin_runs_moments(p, q, l, n; mode="length_ge")
        @test isapprox(m_lge.mean, exact_len_ge; atol=1e-10)
        @test isapprox(m_lge.var, exact_len_ge_sq - exact_len_ge^2; atol=1e-10)

        m_ceq = weighted_bin_runs_moments(p, q, l, n; mode="count_eq")
        @test isapprox(m_ceq.mean, exact_count_eq; atol=1e-10)
        @test isapprox(m_ceq.var, exact_count_eq_sq - exact_count_eq^2; atol=1e-10)

        m_dge = weighted_bin_runs_moments(p, q, l, n; mode="double_count_ge")
        @test isapprox(m_dge.mean, exact_double_ge; atol=1e-10)
        @test isapprox(m_dge.var, exact_double_ge_sq - exact_double_ge^2; atol=1e-10)
    end

    @testset "Markov Binary Runs (Phase 1.1)" begin
        p00, p01, p10, p11 = 0.7, 0.3, 0.2, 0.8
        p_stat = p01 / (p01 + p10)
        q_stat = 1.0 - p_stat
        l, n = 3, 7

        # Brute-force enumeration over all 2^n binary words with 1st-order Markov weights
        ex_cge = 0.0; ex_cge2 = 0.0
        ex_lge = 0.0; ex_lge2 = 0.0
        ex_ceq = 0.0; ex_ceq2 = 0.0
        ex_dge = 0.0; ex_dge2 = 0.0
        ex_no_run = 0.0

        for mask in 0:(2^n - 1)
            bits = [(mask >> k) & 1 for k in 0:(n - 1)]
            w = bits[1] == 1 ? p_stat : q_stat
            for idx in 1:(n - 1)
                a, b = bits[idx], bits[idx + 1]
                w *= (a == 0 ? (b == 0 ? p00 : p01) : (b == 0 ? p10 : p11))
            end

            runs_1 = Int[]
            runs_all = Int[]
            cur_val = bits[1]
            cur_len = 1
            for idx in 2:n
                if bits[idx] == cur_val
                    cur_len += 1
                else
                    push!(runs_all, cur_len)
                    cur_val == 1 && push!(runs_1, cur_len)
                    cur_val = bits[idx]
                    cur_len = 1
                end
            end
            push!(runs_all, cur_len)
            cur_val == 1 && push!(runs_1, cur_len)

            c_ge = count(r -> r >= l, runs_1)
            l_ge = sum(r for r in runs_1 if r >= l; init=0)
            c_eq = count(r -> r == l, runs_1)
            d_ge = count(r -> r >= l, runs_all)

            ex_cge += w * c_ge; ex_cge2 += w * c_ge^2
            ex_lge += w * l_ge; ex_lge2 += w * l_ge^2
            ex_ceq += w * c_eq; ex_ceq2 += w * c_eq^2
            ex_dge += w * d_ge; ex_dge2 += w * d_ge^2
            c_ge == 0 && (ex_no_run += w)
        end

        mc = markov_bin_runs_moments(p00, p01, p10, p11, l, n; mode="count_ge")
        @test isapprox(mc.mean, ex_cge; atol=1e-10)
        @test isapprox(mc.var, ex_cge2 - ex_cge^2; atol=1e-10)

        ml = markov_bin_runs_moments(p00, p01, p10, p11, l, n; mode="length_ge")
        @test isapprox(ml.mean, ex_lge; atol=1e-10)
        @test isapprox(ml.var, ex_lge2 - ex_lge^2; atol=1e-10)

        me = markov_bin_runs_moments(p00, p01, p10, p11, l, n; mode="count_eq")
        @test isapprox(me.mean, ex_ceq; atol=1e-10)
        @test isapprox(me.var, ex_ceq2 - ex_ceq^2; atol=1e-10)

        md = markov_bin_runs_moments(p00, p01, p10, p11, l, n; mode="double_count_ge")
        @test isapprox(md.mean, ex_dge; atol=1e-10)
        @test isapprox(md.var, ex_dge2 - ex_dge^2; atol=1e-10)

        @test isapprox(markov_bin_runs_prob(p00, p01, p10, p11, l, n), ex_no_run; atol=1e-10)
        @test isapprox(markov_bin_runs_pval(p00, p01, p10, p11, l, n), 1.0 - ex_no_run; atol=1e-10)

        # Reduction to i.i.d. Bernoulli when p01 = p11 = p
        p_iid, q_iid = 0.35, 0.65
        m_mkv_iid = markov_bin_runs_moments(q_iid, p_iid, q_iid, p_iid, 3, 12; mode="length_ge")
        m_ber_iid = weighted_bin_runs_moments(p_iid, q_iid, 3, 12; mode="length_ge")
        @test isapprox(m_mkv_iid.mean, m_ber_iid.mean; atol=1e-10)
        @test isapprox(m_mkv_iid.var, m_ber_iid.var; atol=1e-10)
        @test isapprox(markov_bin_runs_prob(q_iid, p_iid, q_iid, p_iid, 3, 12),
                       weighted_bin_runs_prob(p_iid, q_iid, 3, 12); atol=1e-10)

        est = estimate_markov_params([0, 0, 1, 1, 1, 0, 1, 0])
        @test isapprox(est.p00 + est.p01, 1.0; atol=1e-12)
        @test isapprox(est.p10 + est.p11, 1.0; atol=1e-12)
    end

    @testset "Autocorrelation & Pattern BGFs (Phase 1.2)" begin
        # Note III.39 example: pattern ababa, gamma(z) = 1 + pa*pb*z^2 + (pa*pb)^2*z^4
        probs_ab = Dict('a' => 0.3, 'b' => 0.7)
        z_sym = SymPy.symbols("z")
        g_ababa = autocorrelation_poly(['a', 'b', 'a', 'b', 'a'], probs_ab, z_sym)
        @test isapprox(Float64(SymPy.subs(g_ababa, z_sym => 0.5)),
                       1.0 + 0.21 * 0.25 + (0.21^2) * (0.5^4); atol=1e-12)

        # Brute-force over {1, 2, 3}^6 for self-overlapping pattern [1, 2, 1]
        probs3 = [0.2, 0.5, 0.3]
        pat = [1, 2, 1]
        n = 6
        ex_ov1 = 0.0; ex_ov2 = 0.0
        ex_nov1 = 0.0; ex_nov2 = 0.0
        p_ge = zeros(Float64, 4) # P(X >= 0, 1, 2, 3)

        function _enum_words(len, r)
            words = Vector{Vector{Int}}()
            cur = ones(Int, len)
            function rec(pos)
                if pos > len
                    push!(words, copy(cur))
                    return
                end
                for sym in 1:r
                    cur[pos] = sym
                    rec(pos + 1)
                end
            end
            rec(1)
            return words
        end

        for w_seq in _enum_words(n, 3)
            wt = prod(probs3[s] for s in w_seq)
            c_ov = 0
            for i in 1:(n - length(pat) + 1)
                w_seq[i:(i + length(pat) - 1)] == pat && (c_ov += 1)
            end
            c_nov = 0
            idx = 1
            while idx <= n - length(pat) + 1
                if w_seq[idx:(idx + length(pat) - 1)] == pat
                    c_nov += 1
                    idx += length(pat)
                else
                    idx += 1
                end
            end
            ex_ov1 += wt * c_ov; ex_ov2 += wt * c_ov^2
            ex_nov1 += wt * c_nov; ex_nov2 += wt * c_nov^2
            for k_th in 0:3
                c_ov >= k_th && (p_ge[k_th + 1] += wt)
            end
        end

        m_ov = pattern_occurrences_moments(pat, probs3, n; overlapping=true)
        @test isapprox(m_ov.mean, ex_ov1; atol=1e-10)
        @test isapprox(m_ov.var, ex_ov2 - ex_ov1^2; atol=1e-10)

        m_nov = pattern_occurrences_moments(pat, probs3, n; overlapping=false)
        @test isapprox(m_nov.mean, ex_nov1; atol=1e-10)
        @test isapprox(m_nov.var, ex_nov2 - ex_nov1^2; atol=1e-10)

        @test isapprox(pattern_waiting_prob(pat, probs3, n), 1.0 - p_ge[2]; atol=1e-10)
        @test isapprox(pattern_pval(pat, probs3, 1, n), p_ge[2]; atol=1e-10)
        @test isapprox(pattern_pval(pat, probs3, 2, n), p_ge[3]; atol=1e-10)
        @test isapprox(pattern_pval(pat, probs3, 3, n), p_ge[4]; atol=1e-10)
    end

    @testset "Transfer Matrices, Graphs & Random Mappings (Phase 1.3)" begin
        # 2-state Markov chain marking 1->1 transitions (both short n=8 and binary-exp n=100)
        T2 = [0.7 0.3; 0.25 0.75]
        M2 = [0.0 0.0; 0.0 1.0]
        m_tm8 = transfer_matrix_moments(T2, M2, 8)
        @test isapprox(m_tm8.count, 1.0; atol=1e-10)
        # Expected 1->1 transitions per step is pi_2 * T_{2,2}
        pi2 = 0.3 / (0.3 + 0.25)
        @test isapprox(m_tm8.mean, 8.0 * pi2 * 0.75; atol=1e-10)

        m_tm100 = transfer_matrix_moments(T2, M2, 100)
        @test isapprox(m_tm100.count, 1.0; atol=1e-10)
        @test isapprox(m_tm100.mean, 100.0 * pi2 * 0.75; atol=1e-9)
        @test m_tm100.var > 0.0

        # Symbolic & Taylor1 transfer_matrix_gf
        z_s, u_s = SymPy.symbols("z u")
        gf_sym = transfer_matrix_gf(T2, z_s; mark_matrix=M2, u=u_s)
        @test isapprox(Float64(SymPy.subs(SymPy.subs(gf_sym, z_s => 0.2), u_s => 1.0)), 1.0 / (1.0 - 0.2); atol=1e-10)
        z_t = AnalyticComb.TaylorSeries.Taylor1(Float64, 5)
        gf_tay = transfer_matrix_gf(T2, z_t; mark_matrix=M2, u=1.0)
        @test all(k -> isapprox(AnalyticComb.TaylorSeries.getcoeff(gf_tay, k), 1.0; atol=1e-10), 0:5)

        # Golden-mean shift adjacency matrix A = [1 1; 1 0]: det(I - z A) = 1 - z - z^2
        A_fib = [1 1; 1 0]
        @test graph_walks_count(A_fib, 4; source=1, target=1) == 5
        @test graph_walks_count(A_fib, 3) == 8
        @test graph_walks_count(A_fib, 4) == 13
        cw = graph_closed_walks_gf(A_fib, 5)
        @test cw.closed_walks == [1.0, 3.0, 4.0, 7.0, 11.0] # Lucas numbers
        @test isapprox(cw.zeta_inv_coeffs, [1.0, -1.0, -1.0]; atol=1e-10)
        @test isapprox(cw.spectral_radius, (1 + sqrt(5)) / 2; atol=1e-10)

        # Random mappings: brute force n=3 (27 mappings) and n=4 (256 mappings)
        rm3 = random_mapping_moments(3)
        @test isapprox(rm3.expected_cyclic_nodes, 17.0 / 9.0; atol=1e-12)
        rm2 = random_mapping_moments(2)
        @test isapprox(rm2.expected_cyclic_nodes, 1.5; atol=1e-12)
        @test isapprox(rm2.expected_components, 1.25; atol=1e-12)
        @test isapprox(rm2.prob_connected, 0.75; atol=1e-12)
        # Flajolet & Sedgewick p. 463 (n=100: ~2.97 components, ~12.20 cyclic nodes)
        rm100 = random_mapping_moments(100)
        @test isapprox(rm100.expected_components, 2.97; atol=0.05)
        @test isapprox(rm100.expected_cyclic_nodes, 12.20; atol=0.05)
    end

    @testset "Dominant Singularities (Phase 1.4)" begin
        z0_fib = dominant_pole(z -> 1.0 - z - z^2)
        @test isapprox(z0_fib, (sqrt(5) - 1.0) / 2.0; atol=1e-11)

        rp1 = runs_dominant_pole(0.5, 0.5, 4)
        rp2 = runs_dominant_pole(0.8, 0.2, 4)
        rp3 = runs_dominant_pole(0.95, 0.05, 4)
        # Higher run probability p => faster decay of run-avoidance (larger z0, smaller growth_rate)
        # and smaller coherent persistence pole gap (1/p - 1 -> 0)
        @test 1.0 < rp1.z0 < rp2.z0 < rp3.z0
        @test rp1.growth_rate > rp2.growth_rate > rp3.growth_rate > 0.0
        @test rp1.persistence_gap > rp2.persistence_gap > rp3.persistence_gap > 0.0

        # Markov pole matches Bernoulli pole when p01 = p11 = p
        mp_iid = markov_dominant_pole(0.2, 0.8, 0.2, 0.8, 4)
        @test isapprox(mp_iid.z0, rp2.z0; atol=1e-11)

        # Increasing Markov persistence p11 shrinks persistence_gap toward 0 and increases run-avoidance z0
        mp_low = markov_dominant_pole(0.7, 0.3, 0.6, 0.4, 4)
        mp_high = markov_dominant_pole(0.7, 0.3, 0.05, 0.95, 4)
        @test 1.0 < mp_low.z0 < mp_high.z0
        @test mp_low.persistence_gap > mp_high.persistence_gap > 0.0
    end
    




end
