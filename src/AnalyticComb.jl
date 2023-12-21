module AnalyticComb

#using SymPy
using Reexport
using TaylorSeries

@reexport using SymPy

export
	#basic operators
	SEQ, MSET, PSET, CYC, 
	
	#misc
	stirling_factorial_asym, stirling_catalan_asym, catalan_factorial, general_trees, triangulations, 

	#compositions and partitions
	I_gf, partitions_gf, partitions_asym, primes_composition_asym, 
	restricted_sum_comp_gf, restricted_sum_comp, restricted_sum_part_gf, restricted_sum_part,
	partitions_max_r, partitions_asym, fixed_size_comps, fixed_size_comps_asym,

	#words and regular languages
	words_without_k_run, longest_run_binary_asym, 
	bin_words_with_k_occurences, bin_words_with_k_occurences_constr,
	p_binary_words_doub_runl

include("operators.jl")
include("misc.jl")
include("parts_comps.jl")
include("words_languages.jl")


end
