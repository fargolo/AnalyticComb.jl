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
    



end
