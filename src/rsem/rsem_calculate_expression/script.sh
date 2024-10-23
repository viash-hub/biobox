#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

if [ "$par_strandedness" == 'forward' ]; then
    strandedness='--strandedness forward'
elif [ "$par_strandedness" == 'reverse' ]; then
    strandedness="--strandedness reverse"
else
    strandedness=''
fi

IFS=";" read -ra input <<< $par_input

INDEX=$(find -L $par_index -name "*.grp" | sed 's/\.grp$//')
echo "$INDEX"
unset_if_false=( par_paired par_quiet par_no_bam_output par_sampling_for_bam par_no_qualities 
                 par_alignments par_bowtie2 par_star par_hisat2_hca par_append_names 
                 par_single_cell_prior par_calc_pme par_calc_ci par_phred64_quals 
                 par_solexa_quals par_star_gzipped_read_file par_star_bzipped_read_file 
                 par_star_output_genome_bam par_estimate_rspd par_keep_intermediate_files 
                 par_time par_run_pRSEM par_cap_stacked_chipseq_reads par_sort_bam_by_read_name par_sort_bam_by_coordinate )

for par in ${unset_if_false[@]}; do
    test_val="${!par}"
    [[ "$test_val" == "false" ]] && unset $par
done

rsem-calculate-expression \
    ${par_quiet:+-q} \
    ${par_no_bam_output:+--no-bam-output} \
    ${par_sampling_for_bam:+--sampling-for-bam} \
    ${par_no_qualities:+--no-qualities} \
    ${par_alignments:+--alignments} \
    ${par_bowtie2:+--bowtie2} \
    ${par_star:+--star} \
    ${par_hisat2_hca:+--hisat2-hca} \
    ${par_append_names:+--append-names} \
    ${par_single_cell_prior:+--single-cell-prior} \
    ${par_calc_pme:+--calc-pme} \
    ${par_calc_ci:+--calc-ci} \
    ${par_phred64_quals:+--phred64-quals} \
    ${par_solexa_quals:+--solexa-quals} \
    ${par_star_gzipped_read_file:+--star-gzipped-read-file} \
    ${par_star_bzipped_read_file:+--star-bzipped-read-file} \
    ${par_star_output_genome_bam:+--star-output-genome-bam} \
    ${par_estimate_rspd:+--estimate-rspd} \
    ${par_keep_intermediate_files:+--keep-intermediate-files} \
    ${par_time:+--time} \
    ${par_run_pRSEM:+--run-pRSEM} \
    ${par_cap_stacked_chipseq_reads:+--cap-stacked-chipseq-reads} \
    ${par_sort_bam_by_read_name:+--sort-bam-by-read-name} \
    ${par_sort_bam_by_coordinate:+--sort-bam-by-coordinate} \
    ${par_fai:+--fai "$par_fai"} \
    ${par_seed:+--seed "$par_seed"} \
    ${par_seed_length:+--seed-length "$par_seed_length"} \
    ${par_bowtie_n:+--bowtie-n "$par_bowtie_n"} \
    ${par_bowtie_e:+--bowtie-e "$par_bowtie_e"} \
    ${par_bowtie_m:+--bowtie-m "$par_bowtie_m"} \
    ${par_bowtie_chunkmbs:+--bowtie-chunkmbs "$par_bowtie_chunkmbs"} \
    ${par_bowtie2_mismatch_rate:+--bowtie2-mismatch-rate "$par_bowtie2_mismatch_rate"} \
    ${par_bowtie2_k:+--bowtie2-k "$par_bowtie2_k"} \
    ${par_bowtie2_sensitivity_level:+--bowtie2-sensitivity-level "$par_bowtie2_sensitivity_level"} \
    ${par_tag:+--tag "$par_tag"} \
    ${par_fragment_length_min:+--fragment-length-min "$par_fragment_length_min"} \
    ${par_fragment_length_max:+--fragment-length-max "$par_fragment_length_max"} \
    ${par_fragment_length_mean:+--fragment-length-mean "$par_fragment_length_mean"} \
    ${par_fragment_length_sd:+--fragment-length-sd "$par_fragment_length_sd"} \
    ${par_num_rspd_bins:+--num-rspd-bins "$par_num_rspd_bins"} \
    ${par_gibbs_burnin:+--gibbs-burnin "$par_gibbs_burnin"} \
    ${par_gibbs_number_of_samples:+--gibbs-number-of-samples "$par_gibbs_number_of_samples"} \
    ${par_gibbs_sampling_gap:+--gibbs-sampling-gap "$par_gibbs_sampling_gap"} \
    ${par_ci_credibility_level:+--ci-credibility-level "$par_ci_credibility_level"} \
    ${par_ci_number_of_samples_per_count_vector:+--ci-number-of-samples-per-count-vector "$par_ci_number_of_samples_per_count_vector"} \
    ${par_temporary_folder:+--temporary-folder "$par_temporary_folder"} \
    ${par_chipseq_peak_file:+--chipseq-peak-file "$par_chipseq_peak_file"} \
    ${par_chipseq_target_read_files:+--chipseq-target-read-files "$par_chipseq_target_read_files"} \
    ${par_chipseq_control_read_files:+--chipseq-control-read-files "$par_chipseq_control_read_files"} \
    ${par_chipseq_read_files_multi_targets:+--chipseq-read-files-multi-targets "$par_chipseq_read_files_multi_targets"} \
    ${par_chipseq_bed_files_multi_targets:+--chipseq-bed-files-multi-targets "$par_chipseq_bed_files_multi_targets"} \
    ${par_n_max_stacked_chipseq_reads:+--n-max-stacked-chipseq-reads "$par_n_max_stacked_chipseq_reads"} \
    ${par_partition_model:+--partition-model "$par_partition_model"} \
    $strandedness \
    ${par_paired:+--paired-end} \
    ${input[*]} \
    $INDEX \
    $par_id
   
[[ -f "${par_id}.genes.results" ]] && mv "${par_id}.genes.results" $par_counts_gene
[[ -f "${par_id}.isoforms.results" ]] && mv "${par_id}.isoforms.results" $par_counts_transcripts
[[ -d "${par_id}.stat" ]] && mv "${par_id}.stat" $par_stat
[[ -f "${par_id}.log" ]] && mv "${par_id}.log" $par_logs
[[ -f "${par_id}.STAR.genome.bam" ]] && mv "${par_id}.STAR.genome.bam" $par_bam_star
[[ -f "${par_id}.genome.bam" ]] && mv "${par_id}.genome.bam" $par_bam_genome
[[ -f "${par_id}.transcript.bam" ]] && mv "${par_id}.transcript.bam" $par_bam_transcript