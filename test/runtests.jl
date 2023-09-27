using AnalyticComb
using Test

@testset "AnalyticComb.jl" begin

    results = map(i->p_binary_words_runl(i,200),3:12) 
    
    # Check if any of the results presents error > 0.01
    @test sum(results - map(x-> round(x;digits=3),[6.549422864002916e-8,
    0.0007075224314585063,
    0.033979610958548824,
    0.1660282049825555,
    0.2574174128392719,
    0.22352195216110582,
    0.14594875198080856,
    0.08292457405578264,
    0.0440234813776148,
    0.022607888203734293]) .< 0.01) != 0

    @test restricted_sum_part(99,[1,5,10,25]) == 213 # Polya problem for partition with restricted sum

end
