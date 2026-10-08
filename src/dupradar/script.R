library("dupRadar")

## VIASH START
par <- list(
  input_bam = "test_data/sample.bam",
  input_gtf = "test_data/genes.gtf",
  id = "test",
  strandedness = 0L,
  paired = FALSE,
  output_dupmatrix = "dup_matrix.txt",
  output_dup_intercept_mqc = "dup_intercept_mqc.txt",
  output_duprate_exp_boxplot = "duprate_exp_boxplot.pdf",
  output_duprate_exp_densplot = "duprate_exp_densityplot.pdf",
  output_duprate_exp_denscurve_mqc = "duprate_exp_density_curve_mqc.txt",
  output_expression_histogram = "expression_hist.pdf",
  output_intercept_slope = "intercept_slope.txt"
)
meta <- list(
  cpus = 1L
)
## VIASH END

input_bam <- par$input_bam
input_gtf <- par$input_gtf
id <- par$id
stranded <- par$strandedness
paired <- par$paired
threads <- meta$cpus %||% 1L

# Fall back to the original dupRadar file names when no path is given
output_dupmatrix <- par$output_dupmatrix %||%
  paste0(id, "_dupMatrix.txt")
output_intercept_mqc <- par$output_dup_intercept_mqc %||%
  paste0(id, "_dup_intercept_mqc.txt")
output_boxplot <- par$output_duprate_exp_boxplot %||%
  paste0(id, "_duprateExpBoxplot.pdf")
output_densplot <- par$output_duprate_exp_densplot %||%
  paste0(id, "_duprate_exp_densplot.pdf")
output_denscurve_mqc <- par$output_duprate_exp_denscurve_mqc %||%
  paste0(id, "_duprateExpDensCurve_mqc.txt")
output_histogram <- par$output_expression_histogram %||%
  paste0(id, "_expressionHist.pdf")
output_intercept_slope <- par$output_intercept_slope %||%
  paste0(id, "_intercept_slope.txt")

if (!grepl("\\.bam$", input_bam)) {
  stop("--input_bam must be a BAM file (.bam), got: ", input_bam)
}

# Log parameters (stderr)
strand_names <- c("unstranded", "forward", "reverse")
message("> Input BAM:     ", input_bam)
message("> Input GTF:     ", input_gtf)
message("> Sample ID:     ", id)
message("> Strandedness:  ", strand_names[stranded + 1])
message("> Paired-end:    ", paired)
message("> Threads:       ", threads)

# Duplicate stats
message(">> Running analyzeDuprates")
dm <- analyzeDuprates(input_bam, input_gtf, stranded, paired, threads)

# dupRadar only reads the duplicate flag (0x400); it does not set it.
# If no read in any gene is flagged, every duplication rate is 0 and the
# plots are meaningless. Mark duplicates upstream (gatk4_markduplicates).
if (all(dm$dupRate %in% c(0, NA))) {
  warning(
    "No duplicate-flagged reads found in ", input_bam, ". ",
    "dupRadar expects a duplicate-marked BAM; all duplication rates are 0. ",
    "Mark duplicates first, e.g. with the biobox component ",
    "gatk4/gatk4_markduplicates.",
    call. = FALSE
  )
}
message(">> Writing duplication matrix to ", output_dupmatrix)
write.table(
  dm,
  file = output_dupmatrix, quote = FALSE, row.names = FALSE, sep = "\t"
)

# Open a pdf device, draw a dupRadar plot, and always close the device
save_plot <- function(path, plot_fun, main) {
  message(">> Writing ", main, " to ", path)
  pdf(path)
  on.exit(dev.off())
  plot_fun(DupMat = dm)
  title(main)
  mtext(id, side = 3)
}

# 2D density scatter plot
save_plot(output_densplot, duprateExpDensPlot, "Density scatter plot")
# Distribution of expression box plot
save_plot(
  output_boxplot, duprateExpBoxplot, "Percent Duplication by Expression"
)
# Distribution of RPK values per gene
save_plot(
  output_histogram, expressionHist, "Distribution of RPK values per gene"
)

message(">> Running duprateExpFit")
fit <- duprateExpFit(DupMat = dm)
message(">> Writing intercept and slope to ", output_intercept_slope)
cat(
  paste("- dupRadar Int (duprate at low read counts):", fit$intercept),
  paste("- dupRadar Sl (progression of the duplication rate):", fit$slope),
  fill = TRUE, labels = id,
  file = output_intercept_slope, append = FALSE
)

# Create a multiqc file dupInt
sample_name <- gsub("Aligned.sortedByCoord.out.markDups", "", id)
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

message(">> Writing MultiQC intercept table to ", output_intercept_mqc)
write(line, file = output_intercept_mqc, append = FALSE)
write(
  paste(sample_name, fit$intercept),
  file = output_intercept_mqc, append = TRUE
)

# Get numbers from dupRadar GLM
curve_x <- sort(log10(dm$RPK))
curve_y <- 100 * predict(fit$glm, data.frame(x = curve_x), type = "response")
# Remove all of the infinite values
keep <- is.finite(curve_x)
curve_x <- curve_x[keep]
curve_y <- curve_y[keep]
# Reduce number of data points
curve_x <- curve_x[seq(1, length(curve_x), 10)]
curve_y <- curve_y[seq(1, length(curve_y), 10)]
# Convert x values back to real counts
curve_x <- 10^curve_x
# Write to file
# The MultiQC header must stay verbatim, so its long lines are not wrapped
# nolint start: line_length_linter.
line <- "#id: dupradar
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
# nolint end

message(">> Writing MultiQC density curve to ", output_denscurve_mqc)
write(line, file = output_denscurve_mqc, append = FALSE)
write.table(
  cbind(curve_x, curve_y),
  file = output_denscurve_mqc,
  quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE
)
message("> Done running dupRadar")
