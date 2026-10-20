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
    



end
