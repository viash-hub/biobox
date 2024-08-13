#!/bin/bash

set -e

## VIASH START
## VIASH END
unset_if_false=(
    par_discard_orphans
    par_ont
    par_seq_bias
    par_gc_bias
    par_pos_bias
    par_meta
    par_discard_orphans_quasi
    par_disable_chaining_heuristic
    par_allow_dovetail
    par_recover_orphans
    par_mimicBT2
    par_mimic_strictBT2
    par_softclip
    par_softclip_overhangs
    par_full_length_alignment
    par_hard_filter
    par_write_mappings
    par_write_qualities
    par_alternative_init_mode
    par_skip_quant
    par_dump_eq
    par_dump_eq_weights
    par_reduce_GC_memory
    par_init_uniform
    par_no_length_correction
    par_no_effective_length_correction
    par_no_single_frag_prob
    par_no_frag_length_dist
    par_no_bias_length_threshold
    par_useEM
    par_useVBOpt
    par_no_Gamma_draw
    par_bootstrap_reproject
    par_quiet
    par_per_transcript_prior
    par_per_nucleotide_prior
    par_write_orphan_links
    par_write_unmapped_names
    par_no_error_model
    par_sample_out
    par_sample_unaligned
    par_gencode
)

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

IFS=";" read -ra unmated_reads <<< $par_unmated_reads
IFS=";" read -ra mates1 <<< $par_mates1
IFS=";" read -ra mates2 <<< $par_mates2
IFS=";" read -ra alignment <<< $par_alignments

salmon quant \
    ${par_lib_type:+-l "${par_lib_type}"} \
    ${par_index:+-i "${par_index}"} \
    ${par_unmated_reads:+-r ${unmated_reads[*]}} \
    ${par_mates1:+-1 ${mates1[*]}} \
    ${par_mates2:+-2 ${mates2[*]}} \
    ${par_alignments:+-a ${alignment[*]}} \
    ${par_discard_orphans:+--discardOrphans} \
    ${par_eqclasses:+-e "${par_eqclasses}"} \
    ${par_targets:+-t "${par_targets}"} \
    ${par_ont:+--ont} \
    ${par_output:+-o "${par_output}"} \
    ${par_seq_bias:+--seqBias} \
    ${par_gc_bias:+--gcBias} \
    ${par_pos_bias:+--posBias} \
    ${meta_cpus:+-p "${meta_cpus}"} \
    ${par_incompat_prior:+--incompatPrior "${par_incompat_prior}"} \
    ${par_gene_map:+-g "${par_gene_map}"} \
    ${par_aux_target_file:+--auxTargetFile "${par_aux_target_file}"} \
    ${par_meta:+--meta} \
    ${par_score_exp:+--scoreExp "${par_score_exp}"} \
    ${par_discard_orphans_quasi:+--discardOrphansQuasi} \
    ${par_consensus_slack:+--consensusSlack "${par_consensus_slack}"} \
    ${par_pre_merge_chain_sub_thresh:+--preMergeChainSubThresh "${par_pre_merge_chain_sub_thresh}"} \
    ${par_post_merge_chain_sub_thresh:+--postMergeChainSubThresh "${par_post_merge_chain_sub_thresh}"} \
    ${par_orphan_chain_sub_thresh:+--orphanChainSubThresh "${par_orphan_chain_sub_thresh}"} \
    ${par_min_score_fraction:+--minScoreFraction "${par_min_score_fraction}"} \
    ${par_mismatch_seed_skip:+--mismatchSeedSkip "${par_mismatch_seed_skip}"} \
    ${par_disable_chaining_heuristic:+--disableChainingHeuristic} \
    ${par_decoy_threshold:+--decoyThreshold "${par_decoy_threshold}"} \
    ${par_ma:+--ma "${par_ma}"} \
    ${par_mp:+--mp "${par_mp}"} \
    ${par_go:+--go "${par_go}"} \
    ${par_ge:+--ge "${par_ge}"} \
    ${par_bandwidth:+--bandwidth "${par_bandwidth}"} \
    ${par_allow_dovetail:+--allowDovetail} \
    ${par_recover_orphans:+--recoverOrphans} \
    ${par_mimicBT2:+--mimicBT2} \
    ${par_mimic_strictBT2:+--mimicStrictBT2} \
    ${par_softclip:+--softclip} \
    ${par_softclip_overhangs:+--softclipOverhangs} \
    ${par_full_length_alignment:+--fullLengthAlignment} \
    ${par_hard_filter:+--hardFilter} \
    ${par_min_aln_prob:+--minAlnProb "${par_min_aln_prob}"} \
    ${par_write_mappings:+--write_mappings="${par_mappings_sam}"} \
    ${par_write_qualities:+--writeQualities} \
    ${par_hit_filter_policy:+--hitFilterPolicy "${par_hit_filter_policy}"} \
    ${par_alternative_init_mode:+--alternativeInitMode} \
    ${par_aux_dir:+--auxDir "${par_aux_dir}"} \
    ${par_skip_quant:+--skipQuant} \
    ${par_dump_eq:+--dumpEq} \
    ${par_dump_eq_weights:+-d "${par_dump_eq_weights}"} \
    ${par_min_assigned_frags:+--minAssignedFrags "${par_min_assigned_frags}"} \
    ${par_reduce_GC_memory:+--reduceGCMemory} \
    ${par_bias_speed_samp:+--biasSpeedSamp "${par_bias_speed_samp}"} \
    ${par_fld_max:+--fldMax "${par_fld_max}"} \
    ${par_fld_mean:+--fldMean "${par_fld_mean}"} \
    ${par_fld_SD:+--fldSD "${par_fld_SD}"} \
    ${par_forgetting_factor:+-f "${par_forgetting_factor}"} \
    ${par_init_uniform:+--initUniform} \
    ${par_max_occs_per_hit:+--maxOccsPerHit "${par_max_occs_per_hit}"} \
    ${par_max_read_occ:+-w "${par_max_read_occ}"} \
    ${par_no_length_correction:+--noLengthCorrection} \
    ${par_no_effective_length_correction:+--noEffectiveLengthCorrection} \
    ${par_no_single_frag_prob:+--noSingleFragProb} \
    ${par_no_frag_length_dist:+--noFragLengthDist} \
    ${par_no_bias_length_threshold:+--noBiasLengthThreshold} \
    ${par_num_bias_samples:+--numBiasSamples "${par_num_bias_samples}"} \
    ${par_num_aux_model_samples:+--numAuxModelSamples "${par_num_aux_model_samples}"} \
    ${par_num_pre_aux_model_samples:+--numPreAuxModelSamples "${par_num_pre_aux_model_samples}"} \
    ${par_useEM:+--useEM} \
    ${par_useVBOpt:+--useVBOpt} \
    ${par_range_factorization_bins:+--rangeFactorizationBins "${par_range_factorization_bins}"} \
    ${par_num_Gibbs_samples:+--numGibbsSamples "${par_num_Gibbs_samples}"} \
    ${par_no_Gamma_draw:+--noGammaDraw} \
    ${par_num_bootstraps:+--numBootstraps "${par_num_bootstraps}"} \
    ${par_bootstrap_reproject:+--bootstrapReproject} \
    ${par_thinning_factor:+--thinningFactor "${par_thinning_factor}"} \
    ${par_quiet:+--quiet} \
    ${par_per_transcript_prior:+--perTranscriptPrior} \
    ${par_per_nucleotide_prior:+--perNucleotidePrior} \
    ${par_sig_digits:+--sigDigits "${par_sig_digits}"} \
    ${par_vb_prior:+--vbPrior "${par_vb_prior}"} \
    ${par_write_orphan_links:+--writeOrphanLinks} \
    ${par_write_unmapped_names:+--writeUnmappedNames} \
    ${par_no_error_model:+--noErrorModel} \
    ${par_num_error_bins:+--numErrorBins "${par_num_error_bins}"} \
    ${par_sample_out:+--sampleOut} \
    ${par_sample_unaligned:+--sampleUnaligned} \
    ${par_gencode:+--gencode} \
    ${par_mapping_cache_memory_limit:+--mappingCacheMemoryLimit "${par_mapping_cache_memory_limit}"}

if [ -f "$par_output/quant.sf" ]; then
    mv $par_output/quant.sf $par_quant_results
else
    echo "Quantification file not generated!"
fi