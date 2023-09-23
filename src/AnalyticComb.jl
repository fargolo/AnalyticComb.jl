module AnalyticComb

using SymPy
#using TaylorSeries
#using DynamicalSystems

export
	stirling_factorial, partitions, primes_composition
	p_binary_words_runl

include("binary_words.jl")
include("asymptotics.jl")

end
