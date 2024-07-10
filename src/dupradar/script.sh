#!/bin/bash
## VIASH START
## VIASH END


set -eo pipefail 

function num_strandness {
    if [ $par_strandedness == 'unstranded' ]; then echo 0
    elif [ $par_strandedness == 'forward' ]; then echo 1
    elif [ $par_strandedness == 'reverse' ]; then echo 2
    else echo "strandedness must be unstranded, forward or reverse." && \
        exit 1
    fi
}

Rscript "$meta_resources_dir/script.R" \
    $par_input \
    $par_id \
    $par_gtf_annotation \
    $(num_strandness) \
    $par_paired \
    ${meta_cpus:-1}

mv "$par_id"_dupMatrix.txt $par_output_dupmatrix
mv "$par_id"_dup_intercept_mqc.txt $par_output_dup_intercept_mqc
mv "$par_id"_duprateExpBoxplot.pdf $par_output_duprate_exp_boxplot
mv "$par_id"_duprateExpDens.pdf $par_output_duprate_exp_densplot
mv "$par_id"_duprateExpDensCurve_mqc.txt $par_output_duprate_exp_denscurve_mqc
mv "$par_id"_expressionHist.pdf $par_output_expression_histogram
mv "$par_id"_intercept_slope.txt $par_output_intercept_slope


dupradar_ver=$(Rscript -e "library(dupRadar); cat(as.character(packageVersion('dupRadar')))")
text="bioconductor-dupradar: ${dupradar_ver}"
echo "$text"

exit 0
