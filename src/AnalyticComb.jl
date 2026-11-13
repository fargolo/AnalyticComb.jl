module AnalyticComb

#using SymPy
using Reexport
using TaylorSeries
using LinearAlgebra

@reexport using SymPy

export
	#basic operators, BGF moments, and dominant singularities
	SEQ, MSET, PSET, CYC, bgf_moments,
	dominant_pole, runs_dominant_pole, markov_dominant_pole,
	
	#misc
	stirling_factorial_asym, stirling_catalan_asym, catalan_factorial, general_trees, triangulations, 

	#compositions and partitions
	I_gf, partitions_gf, partitions_asym, primes_composition_asym, 
	restricted_sum_comp_gf, restricted_sum_comp, restricted_sum_part_gf, restricted_sum_part,
	partitions_max_r, partitions_asym, fixed_size_comps, fixed_size_comps_asym,

	#words, regular languages, Markov runs & arbitrary patterns
	words_without_k_run, longest_run_binary_asym, 
	bin_words_with_k_occurences, bin_words_with_k_occurences_constr,
	bin_words_runs_prob, weighted_bin_runs_pval , weighted_bin_runs_prob,
	weighted_bin_runs_bgf, weighted_bin_runs_moments,
	estimate_markov_params, markov_bin_runs_bgf, markov_bin_runs_moments,
	markov_bin_runs_prob, markov_bin_runs_pval,
	autocorrelation_poly, pattern_occurrences_bgf, pattern_occurrences_moments,
	pattern_waiting_prob, pattern_pval,

	#transfer matrices, graphs, dynamical zeta & random mappings
	transfer_matrix_moments, transfer_matrix_gf,
	graph_walks_count, graph_closed_walks_gf, random_mapping_moments

include("operators.jl")
include("misc.jl")
include("parts_comps.jl")
include("words_languages.jl")
include("transfer_graphs.jl")


end
