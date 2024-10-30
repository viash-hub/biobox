#!/usr/bin/env Rscript

## VIASH START
par <- list(
  "input" = 'test_data/sample.bam',
  "id" = 'test',
  "gtf_annotation" = 'test_data/genes.gtf',
  "strandedness" = 1,
  "paired" = TRUE,
  "threads" = 1,
  "output_dupmatrix"="dup_matrix.txt",
  "output_dup_intercept_mqc"="dup_intercept_mqc.txt",
  "output_duprate_exp_boxplot"="duprate_exp_boxplot.pdf",
  "output_duprate_exp_densplot"="duprate_exp_densityplot.pdf",
  "output_duprate_exp_denscurve_mqc"="duprate_exp_density_curve_mqc.pdf",
  "output_expression_histogram"="expression_hist.pdf",
  "output_intercept_slope"="intercept_slope.txt"
)
## VIASH END

library(rlang)

if (length(par) < 5) {
    stop("Usage: dupRadar.r <input.bam> <sample_id> <annotation.gtf> <strandDirection:0=unstranded/1=forward/2=reverse> <paired/single> <nbThreads> <R-package-location (optional)>", call.=FALSE)
}


input_bam <- par$input
output_prefix <- par$id
annotation_gtf <- par$gtf_annotation
stranded <- par$strandedness %||% 0L
paired_end <- ifelse(par$paired, TRUE, FALSE)
threads <- ifelse(is.null(par$threads), 1, as.integer(par$threads))

output_dupmatrix <- par$output_dupmatrix %||% paste0(output_prefix, "_dupMatrix.txt")
output_dup_intercept_mqc <- ifelse(is.null(par$output_dup_intercept_mqc), paste0(output_prefix, "_dup_intercept_mqc.txt", sep=""), par$output_dup_intercept_mqc)
output_duprate_exp_boxplot <- ifelse(is.null(par$output_duprate_exp_boxplot), paste0(output_prefix, "_duprateExpBoxplot.pdf", sep=""), par$output_duprate_exp_boxplot)
output_duprate_exp_densplot <- ifelse(is.null(par$output_duprate_exp_densplot), paste0(output_prefix, "_duprate_exp_densplot.pdf", sep=""), par$output_duprate_exp_densplot)
output_duprate_exp_denscurve_mqc <- ifelse(is.null(par$output_duprate_exp_denscurve_mqc), paste0(output_prefix, "_duprateExpDensCurve_mqc.txt", sep=""), par$output_duprate_exp_denscurve_mqc)
output_expression_histogram <- ifelse(is.null(par$output_expression_histogram), paste0(output_prefix, "_expressionHist.pdf", sep=""), par$output_expression_histogram)
output_intercept_slope <- ifelse(is.null(par$output_intercept_slope), paste0(output_prefix, "_intercept_slope.txt", sep=""), par$output_intercept_slope)

bamRegex <- "(.+)\\.bam$"


if(!(grepl(bamRegex, input_bam) && file.exists(input_bam) && (!file.info(input_bam)$isdir))) stop("First argument '<input.bam>' must be an existing file (not a directory) with '.bam' extension...")
if(!(file.exists(annotation_gtf) &&  (!file.info(annotation_gtf)$isdir))) stop("Third argument '<annotation.gtf>' must be an existing file (and not a directory)...")
if(is.na(stranded) || (!(stranded %in% (0:2)))) stop("Fourth argument <strandDirection> must be a numeric value in 0(unstranded)/1(forward)/2(reverse)...")
if(is.na(threads) || (threads<=0)) stop("Fifth argument <nbThreads> must be a strictly positive numeric value...")

# Debug messages (stderr)
message("Input bam      (Arg 1): ", input_bam)
message("Output basename(Arg 2): ", output_prefix)
message("Input gtf      (Arg 3): ", annotation_gtf)
message("Strandness     (Arg 4): ", c("unstranded", "forward", "reverse")[stranded])
message("paired_end     (Arg 5): ", paired_end)
message("Nb threads     (Arg 6): ", threads)
message("R package loc. (Arg 7): ", ifelse(length(args) > 4, args[5], "Not specified"))


# Load / install packages
if (length(args) > 5) { .libPaths( c( args[6], .libPaths() ) ) }
if (!require("dupRadar")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("dupRadar", update = TRUE, ask=FALSE)
  library("dupRadar")
}
if (!require("parallel")) {
    library("parallel")
}


# Duplicate stats
dm <- analyzeDuprates(input_bam, annotation_gtf, stranded, paired_end, threads)
print("analyzeDuprates done")
write.table(dm, file=output_dupmatrix, quote=F, row.name=F, sep="\t")
print("write.table done")

# 2D density scatter plot
pdf(output_duprate_exp_densplot)
print("pdf done")
duprateExpDensPlot(DupMat=dm)
title("Density scatter plot")
mtext(output_prefix, side=3)
dev.off()
print("duprateExpDensPlot done")
fit <- duprateExpFit(DupMat=dm)
print("duprateExpFit done")
cat(
    paste("- dupRadar Int (duprate at low read counts):", fit$intercept),
    paste("- dupRadar Sl (progression of the duplication rate):", fit$slope),
    fill=TRUE, labels=output_prefix,
    file=output_intercept_slope, append=FALSE
)

# Create a multiqc file dupInt
sample_name <- gsub("Aligned.sortedByCoord.out.markDups", "", output_prefix)
line <- "#id: DupInt
#plot_type: 'generalstats'
#pconfig:
#    dupRadar_intercept:
#        title: 'dupInt'
#        namespace: 'DupRadar'
#        description: 'Intercept value from DupRadar'
#        max: 100
#        min: 0
#        scale: 'RdYlGn-rev'
#        format: '{:.2f}%'
Sample dupRadar_intercept"

write(line,file=output_dup_intercept_mqc, append=TRUE)
write(paste(sample_name, fit$intercept),file=output_dup_intercept_mqc, append=TRUE)
print("write dup_intercept_mqc done")

# Get numbers from dupRadar GLM
curve_x <- sort(log10(dm$RPK))
curve_y = 100*predict(fit$glm, data.frame(x=curve_x), type="response")
# Remove all of the infinite values
infs = which(curve_x %in% c(-Inf,Inf))
curve_x = curve_x[-infs]
curve_y = curve_y[-infs]
# Reduce number of data points
curve_x <- curve_x[seq(1, length(curve_x), 10)]
curve_y <- curve_y[seq(1, length(curve_y), 10)]
# Convert x values back to real counts
curve_x = 10^curve_x
# Write to file
line="#id: dupradar
#section_name: 'DupRadar'
#section_href: 'bioconductor.org/packages/release/bioc/html/dupRadar.html'
#description: \"provides duplication rate quality control for RNA-Seq datasets. Highly expressed genes can be expected to have a lot of duplicate reads, but high numbers of duplicates at low read counts can indicate low library complexity with technical duplication.
#    This plot shows the general linear models - a summary of the gene duplication distributions. \"
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xLog: True
#    xlab: 'expression (reads/kbp)'
#    ylab: '% duplicate reads'
#    ymax: 100
#    ymin: 0
#    tt_label: '<b>{point.x:.1f} reads/kbp</b>: {point.y:,.2f}% duplicates'
#    xPlotLines:
#        - color: 'green'
#          dashStyle: 'LongDash'
#          label:
#                style: {color: 'green'}
#                text: '0.5 RPKM'
#                verticalAlign: 'bottom'
#                y: -65
#          value: 0.5
#          width: 1
#        - color: 'red'
#          dashStyle: 'LongDash'
#          label:
#                style: {color: 'red'}
#                text: '1 read/bp'
#                verticalAlign: 'bottom'
#                y: -65
#          value: 1000
#          width: 1"

write(line,file=output_duprate_exp_denscurve_mqc, append=TRUE)
write.table(
    cbind(curve_x, curve_y),
    file=output_duprate_exp_denscurve_mqc,
    quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE,
)
print("write duprateExpDensCurve_mqc done")
# Distribution of expression box plot
pdf(output_duprate_exp_boxplot)
duprateExpBoxplot(DupMat=dm)
title("Percent Duplication by Expression")
mtext(output_prefix, side=3)
dev.off()
print("duprateExpBoxplot done")
# Distribution of RPK values per gene
pdf(output_expression_histogram)
expressionHist(DupMat=dm)
title("Distribution of RPK values per gene")
mtext(output_prefix, side=3)
dev.off()
print("expressionHist done")


# 

# Print sessioninfo to standard out
print(output_prefix)
# citation("dupRadar")
# print("citation done")
# sessionInfo()
