using AnalyticComb
using Test

@testset "AnalyticComb.jl" begin

    @testset "Words" begin
       
        # Binary words with distance constraints
        ## 4 occurences in a 10-letter word separated by maximum distance of 2 
        @test bin_words_with_k_occurences_constr(4,10,2) == 45

        @test words_without_k_run(3,10) == 504

        binary_double_results = map(i->p_binary_words_doub_runl(i,200),3:12)
             
        # Check if any of the results presents error > 0.01
        @test sum(binary_double_results - map(x-> round(x;digits=3),
            [6.5494e-8,0.00070,0.03397,0.16602,0.25741,0.22352,0.14594,0.08292,0.04402,0.02260]) .< 0.01) != 0
        
    end
    
    # Polya problem for partition with restricted sums (partitions) for 99. 
    @test restricted_sum_part(99,[1,5,10,25]) == 213 


end
