import subprocess
import gzip
from pathlib import Path
from typing import Tuple
import numpy as np
import random
import mudata as md

## VIASH START
meta = {
  "name": "bd_rhapsody_sequence_analysis",
  "executable": "target/docker/bd_rhapsody/bd_rhapsody_sequence_analysis/bd_rhapsody_sequence_analysis",
  "resources_dir": "src/bd_rhapsody",
  "cpus": 8,
  "memory_mb": 4096,
}
## VIASH END

import sys
sys.path.append(meta["resources_dir"])

from helpers.rhapsody_cell_label import index_to_sequence

meta["executable"] = Path(meta["executable"])
meta["resources_dir"] = Path(meta["resources_dir"])

#########################################################################################

# Generate index
print("> Generate index", flush=True)
# cwl_file = meta["resources_dir"] / "bd_rhapsody_make_reference.cwl"
cwl_file = "/var/bd_rhapsody_cwl/v2.2.1/Extra_Utilities/make_rhap_reference_2.2.1.cwl"
reference_small_gtf = meta["resources_dir"] / "test_data" / "reference_small.gtf"
reference_small_fa = meta["resources_dir"] / "test_data" / "reference_small.fa"
bdabseq_panel_fa = meta["resources_dir"] / "test_data" / "BDAbSeq_ImmuneDiscoveryPanel.fasta"
sampletagsequences_fa = meta["resources_dir"] / "test_data" / "SampleTagSequences_HomoSapiens_ver1.fasta"

config_file = Path("reference_config.yml")
reference_file = Path("Rhap_reference.tar.gz")

subprocess.run([
    "cwl-runner", 
    "--no-container",
    "--preserve-entire-environment",
    "--outdir",
    ".",
    str(cwl_file),
    "--Genome_fasta",
    str(reference_small_fa),
    "--Gtf", 
    str(reference_small_gtf),
    "--Extra_STAR_params",
    "--genomeSAindexNbases 4"
])

#########################################################################################
# Load reference in memory

from Bio import SeqIO
import gffutils

# Load FASTA sequence
with open(str(reference_small_fa), "r") as handle:
  reference_fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
with open(str(bdabseq_panel_fa), "r") as handle:
  bdabseq_panel_fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
with open(str(sampletagsequences_fa), "r") as handle:
  sampletagsequences_fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

# create in memory db
reference_gtf_db = gffutils.create_db(
  str(reference_small_gtf),
  dbfn=":memory:",
  force=True,
  keep_order=True,
  merge_strategy="merge",
  sort_attribute_values=True,
  disable_infer_transcripts=True,
  disable_infer_genes=True
)

#############################################
# TODO: move helper functions to separate helper file


def generate_bd_read_metadata(
  instrument_id: str = "A00226",
  run_id: str = "970",
  flowcell_id: str = "H5FGVMXY",
  lane: int = 1,
  tile: int = 1101,
  x: int = 1000,
  y: int = 1000,
  illumina_flag: str = "1:N:0",
  sample_id: str = "CAGAGAGG",
) -> str:
  """
  Generate a FASTQ metadata line for a BD Rhapsody FASTQ file.

  Args:
    instrument_id: The instrument ID.
    run_id: The run ID.
    flowcell_id: The flowcell ID.
    lane: The lane number.
    tile: The tile number. Between 1101 and 1112 in the used example data.
    x: The x-coordinate. Between 1000 and 32967 in the used example data.
    y: The y-coordinate. Between 1000 and 37059 in the used example data.
    illumina_flag: The Illumina flag. Either 1:N:0 or 2:N:0 in the used example data.
    sample_id: The sample ID.
  """
  # format: @A00226:970:H5FGVDMXY:1:1101:2645:1000 2:N:0:CAGAGAGG
  return f"@{instrument_id}:{run_id}:{flowcell_id}:{lane}:{tile}:{x}:{y} {illumina_flag}:{sample_id}"


def generate_bd_wta_transcript(
  transcript_length: int = 42,
) -> str:
  """
  Generate a WTA transcript from a given GTF and FASTA file.
  """

  # Randomly select a gene
  gene = random.choice(list(reference_gtf_db.features_of_type("gene")))

  # Find all exons within the gene
  exons = list(reference_gtf_db.children(gene, featuretype="exon", order_by="start"))

  # Calculate total exon length
  total_exon_length = sum(exon.end - exon.start + 1 for exon in exons)

  # If total exon length is less than desired transcript length, use it as is
  max_transcript_length = min(total_exon_length, transcript_length)

  # Build the WTA transcript sequence
  sequence = ""
  for exon in exons:
    exon_seq = str(reference_fasta_dict[exon.seqid].seq[exon.start - 1 : exon.end])  
    sequence += exon_seq

    # Break if desired length is reached
    if len(sequence) >= max_transcript_length:
      sequence = sequence[:max_transcript_length]
      break
  
  # add padding if need be
  if len(sequence) < max_transcript_length:
    sequence += "N" * (max_transcript_length - len(sequence))

  return sequence


def generate_bd_wta_read(
  cell_index: int = 0,
  bead_version: str = "EnhV2",
  umi_length: int = 14,
  transcript_length: int = 42,
) -> Tuple[str, str]:
  """
  Generate a BD Rhapsody WTA read pair for a given cell index.

  Args:
    cell_index: The cell index to generate reads for.
    bead_version: The bead version to use for generating the cell label.
    umi_length: The length of the UMI to generate.
    transcript_length: The length of the transcript to generate

  Returns:
    A tuple of two strings, the first string being the R1 read and the second string being the R2 read.

  More info:
    
    See structure of reads:
    - https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/top_steps.html
    - https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/steps_cell_label.html
    - https://scomix.bd.com/hc/en-us/articles/360057714812-All-FAQ
    R1 is Cell Label + UMI + PolyT -> 60 bp
      actually, CLS1 + "GTGA" + CLS2 + "GACA" + CLS3 + UMI
    R2 is the actual read -> 42 bp

    Example R1
    CLS1       Link CLS2      Link CLS3       UMI
    AAAATCCTGT GTGA AACCAAAGT GACA GATAGAGGAG CGCATGTTTATAAC
  """
  
  # generate metadata
  per_row = np.floor((32967 - 1000) / 9)
  per_col = np.floor((37059 - 1000) / 9)
  
  assert cell_index >= 0 and cell_index < per_row * per_col, f"cell_index must be between 0 and {per_row} * {per_col}"
  x = 1000 + (cell_index % per_row) * 9
  y = 1000 + (cell_index // per_row) * 9
  instrument_id = "A00226"
  run_id = "970"
  flowcell_id = "H5FGVMXY"
  meta_r1 = generate_bd_read_metadata(instrument_id=instrument_id, run_id=run_id, flowcell_id=flowcell_id, x=x, y=y, illumina_flag="1:N:0")
  meta_r2 = generate_bd_read_metadata(instrument_id=instrument_id, run_id=run_id, flowcell_id=flowcell_id, x=x, y=y, illumina_flag="2:N:0")

  # generate r1 (cls1 + link + cls2 + link + cls3 + umi)
  assert cell_index >= 0 and cell_index < 384 * 384 * 384
  cell_label = index_to_sequence(cell_index + 1, bead_version=bead_version)
  # sample random umi
  umi = "".join(random.choices("ACGT", k=umi_length))
  quality_r1 = "I" * (len(cell_label) + len(umi))
  r1 = f"{meta_r1}\n{cell_label}{umi}\n+\n{quality_r1}\n"

  # generate r2 by extracting sequence from fasta and gtf
  wta_transcript = generate_bd_wta_transcript(transcript_length=transcript_length)
  quality_r2 = "I" * transcript_length
  r2 = f"{meta_r2}\n{wta_transcript}\n+\n{quality_r2}\n"

  return r1, r2

def generate_bd_wta_fastq_files(
  num_cells: int = 100,
  num_reads_per_cell: int = 1000,
) -> Tuple[str, str]:
  """
  Generate BD Rhapsody WTA FASTQ files for a given number of cells and transcripts per cell.

  Args:
    num_cells: The number of cells to generate
    num_reads_per_cell: The number of reads to generate per cell

  Returns:
    A tuple of two strings, the first string being the R1 reads and the second string being the R2 reads.
  """
  r1_reads = ""
  r2_reads = ""
  for cell_index in range(num_cells):
    for _ in range(num_reads_per_cell):
      r1, r2 = generate_bd_wta_read(cell_index)
      r1_reads += r1
      r2_reads += r2

  return r1_reads, r2_reads

def generate_bd_abc_read(
  cell_index: int = 0,
  bead_version: str = "EnhV2",
  umi_length: int = 14,
  transcript_length: int = 72,
) -> Tuple[str, str]:
  """
  Generate a BD Rhapsody ABC read pair for a given cell index.

  Args:
    cell_index: The cell index to generate reads for.
    bead_version: The bead version to use for generating the cell label.
    umi_length: The length of the UMI to generate.
    transcript_length: The length of the transcript to generate

  Returns:
    A tuple of two strings, the first string being the R1 read and the second string being the R2 read.
  """
  # generate metadata
  per_row = np.floor((32967 - 1000) / 9)
  per_col = np.floor((37059 - 1000) / 9)
  
  assert cell_index >= 0 and cell_index < per_row * per_col, f"cell_index must be between 0 and {per_row} * {per_col}"
  x = 1000 + (cell_index % per_row) * 9
  y = 1000 + (cell_index // per_row) * 9
  instrument_id = "A01604"
  run_id = "19"
  flowcell_id = "HMKLYDRXY"
  meta_r1 = generate_bd_read_metadata(instrument_id=instrument_id, run_id=run_id, flowcell_id=flowcell_id, x=x, y=y, illumina_flag="1:N:0")
  meta_r2 = generate_bd_read_metadata(instrument_id=instrument_id, run_id=run_id, flowcell_id=flowcell_id, x=x, y=y, illumina_flag="2:N:0")

  # generate r1 (cls1 + link + cls2 + link + cls3 + umi)
  assert cell_index >= 0 and cell_index < 384 * 384 * 384
  cell_label = index_to_sequence(cell_index + 1, bead_version=bead_version)
  # sample random umi
  umi = "".join(random.choices("ACGT", k=umi_length))
  quality_r1 = "I" * (len(cell_label) + len(umi))
  r1 = f"{meta_r1}\n{cell_label}{umi}\n+\n{quality_r1}\n"

  # generate r2 by sampling sequence from bdabseq_panel_fa
  abseq_seq = str(random.choice(list(bdabseq_panel_fasta_dict.values())).seq)
  abc_suffix = "AAAAAAAAAAAAAAAAAAAAAAA"
  abc_data = abseq_seq[:transcript_length - len(abc_suffix) - 1]
  abc_prefix = "N" + "".join(random.choices("ACGT", k=transcript_length - len(abc_data) - len(abc_suffix) - 1))

  abc_transcript = f"{abc_prefix}{abc_data}{abc_suffix}"

  quality_r2 = "#" + "I" * (len(abc_transcript) - 1)
  r2 = f"{meta_r2}\n{abc_transcript}\n+\n{quality_r2}\n"

  return r1, r2

def generate_bd_abc_fastq_files(
  num_cells: int = 100,
  num_reads_per_cell: int = 1000,
) -> Tuple[str, str]:
  """
  Generate BD Rhapsody ABC FASTQ files for a given number of cells and transcripts per cell.

  Args:
    num_cells: The number of cells to generate
    num_reads_per_cell: The number of reads to generate per cell

  Returns:
    A tuple of two strings, the first string being the R1 reads and the second string being the R2 reads.
  """
  r1_reads = ""
  r2_reads = ""
  for cell_index in range(num_cells):
    for _ in range(num_reads_per_cell):
      r1, r2 = generate_bd_abc_read(cell_index)
      r1_reads += r1
      r2_reads += r2

  return r1_reads, r2_reads

def generate_bd_smk_read(
  cell_index: int = 0,
  bead_version: str = "EnhV2",
  umi_length: int = 14,
  transcript_length: int = 72,
  num_sample_tags: int = 3,
):
  """
  Generate a BD Rhapsody SMK read pair for a given cell index.

  Args:
    cell_index: The cell index to generate reads for.
    bead_version: The bead version to use for generating the cell label.
    umi_length: The length of the UMI to generate.
    transcript_length: The length of the transcript to generate
    num_sample_tags: The number of sample tags to use

  Returns:
    A tuple of two strings, the first string being the R1 read and the second string being the R2 read.
  """
  # generate metadata
  per_row = np.floor((32967 - 1000) / 9)
  per_col = np.floor((37059 - 1000) / 9)
  
  assert cell_index >= 0 and cell_index < per_row * per_col, f"cell_index must be between 0 and {per_row} * {per_col}"
  x = 1000 + (cell_index % per_row) * 9
  y = 1000 + (cell_index // per_row) * 9
  instrument_id = "A00226"
  run_id = "970"
  flowcell_id = "H5FGVDMXY"
  
  meta_r1 = generate_bd_read_metadata(instrument_id=instrument_id, run_id=run_id, flowcell_id=flowcell_id, x=x, y=y, illumina_flag="1:N:0")
  meta_r2 = generate_bd_read_metadata(instrument_id=instrument_id, run_id=run_id, flowcell_id=flowcell_id, x=x, y=y, illumina_flag="2:N:0")

  # generate r1 (cls1 + link + cls2 + link + cls3 + umi)
  assert cell_index >= 0 and cell_index < 384 * 384 * 384
  cell_label = index_to_sequence(cell_index + 1, bead_version=bead_version)
  # sample random umi
  umi = "".join(random.choices("ACGT", k=umi_length))
  quality_r1 = "I" * (len(cell_label) + len(umi))
  r1 = f"{meta_r1}\n{cell_label}{umi}\n+\n{quality_r1}\n"

  # generate r2 by selecting the cell_index %% num_sample_tags sample tags
  sampletag_index = cell_index % num_sample_tags
  sampletag_seq = str(list(sampletagsequences_fasta_dict.values())[sampletag_index].seq)
  smk_data = sampletag_seq[:transcript_length]
  smk_suffix = "A" * (transcript_length - len(smk_data))
  quality_r2 = "I" * len(smk_data) + "#" * len(smk_suffix)
  r2 = f"{meta_r2}\n{smk_data}{smk_suffix}\n+\n{quality_r2}\n"

  return r1, r2

def generate_bd_smk_fastq_files(
  num_cells: int = 100,
  num_reads_per_cell: int = 1000,
  num_sample_tags: int = 3,
) -> Tuple[str, str]:
  """
  Generate BD Rhapsody SMK FASTQ files for a given number of cells and transcripts per cell.

  Args:
    num_cells: The number of cells to generate
    num_reads_per_cell: The number of reads to generate per cell
    num_sample_tags: The number of sample tags to use

  Returns:
    A tuple of two strings, the first string being the R1 reads and the second string being the R2 reads.
  """
  r1_reads = ""
  r2_reads = ""
  for cell_index in range(num_cells):
    for _ in range(num_reads_per_cell):
      r1, r2 = generate_bd_smk_read(cell_index, num_sample_tags=num_sample_tags)
      r1_reads += r1
      r2_reads += r2

  return r1_reads, r2_reads

#########################################################################################

# Prepare WTA, ABC, and SMK test data
print("> Prepare WTA test data", flush=True)
wta_reads_r1_str, wta_reads_r2_str = generate_bd_wta_fastq_files(num_cells=100, num_reads_per_cell=1000)
with gzip.open("WTAreads_R1.fq.gz", "wt") as f:
  f.write(wta_reads_r1_str)
with gzip.open("WTAreads_R2.fq.gz", "wt") as f:
  f.write(wta_reads_r2_str)

print("> Prepare ABC test data", flush=True)
abc_reads_r1_str, abc_reads_r2_str = generate_bd_abc_fastq_files(num_cells=100, num_reads_per_cell=1000)
with gzip.open("ABCreads_R1.fq.gz", "wt") as f:
  f.write(abc_reads_r1_str)
with gzip.open("ABCreads_R2.fq.gz", "wt") as f:
  f.write(abc_reads_r2_str)

print("> Prepare SMK test data", flush=True)
smk_reads_r1_str, smk_reads_r2_str = generate_bd_smk_fastq_files(num_cells=100, num_reads_per_cell=1000, num_sample_tags=3)
with gzip.open("SMKreads_R1.fq.gz", "wt") as f:
  f.write(smk_reads_r1_str)
with gzip.open("SMKreads_R2.fq.gz", "wt") as f:
  f.write(smk_reads_r2_str)

#########################################################################################

# Run executable
print(f">> Run {meta['name']}", flush=True)
output_dir = Path("output")
subprocess.run([
  meta['executable'],
  "--reads=WTAreads_R1.fq.gz;WTAreads_R2.fq.gz",
  f"--reference_archive={reference_file}",
  "--reads=ABCreads_R1.fq.gz;ABCreads_R2.fq.gz",
  f"--abseq_reference={bdabseq_panel_fa}",
  "--reads=SMKreads_R1.fq.gz;SMKreads_R2.fq.gz",
  "--tag_names=1-Sample1;2-Sample2;3-Sample3",
  "--sample_tags_version=human",
  "--output_dir=output",
  "--exact_cell_count=100",
  f"---cpus={meta['cpus'] or 1}",
  f"---memory={meta['memory_mb'] or 2048}mb",
  # "--output_seurat=seurat.rds",
  "--output_mudata=mudata.h5mu",
  "--metrics_summary=metrics_summary.csv",
  "--pipeline_report=pipeline_report.html",
])


# Check if output exists
print(">> Check if output exists", flush=True)
assert (output_dir / "sample_Bioproduct_Stats.csv").exists()
assert (output_dir / "sample_Metrics_Summary.csv").exists()
assert (output_dir / "sample_Pipeline_Report.html").exists()
assert (output_dir / "sample_RSEC_MolsPerCell_MEX.zip").exists()
assert (output_dir / "sample_RSEC_MolsPerCell_Unfiltered_MEX.zip").exists()
# seurat object is not generated when abc data is added
# assert (output_dir / "sample_Seurat.rds").exists()
assert (output_dir / "sample.h5mu").exists()

# check individual outputs
# assert Path("seurat.rds").exists()
assert Path("mudata.h5mu").exists()
assert Path("metrics_summary.csv").exists()
assert Path("pipeline_report.html").exists()

print(">> Check contents of output", flush=True)
data = md.read_h5mu("mudata.h5mu")

assert data.n_obs == 100, "Number of cells is incorrect"
assert "rna" in data.mod, "RNA data is missing"
assert "prot" in data.mod, "Protein data is missing"

# check rna data
data_rna = data.mod["rna"]
assert data_rna.n_vars == 1, "Number of genes is incorrect"
assert data_rna.X.sum(axis=1).min() > 950, "Number of reads per cell is incorrect"
# assert data_rna.var.Raw_Reads.sum() == 100000, "Number of reads is incorrect"
assert data_rna.var.Raw_Reads.sum() >= 99990 and data_rna.var.Raw_Reads.sum() <= 100010, \
  f"Expected 100000 RNA reads, got {data_rna.var.Raw_Reads.sum()}"

# check prot data
data_prot = data.mod["prot"]
assert data_prot.n_vars == len(bdabseq_panel_fasta_dict), "Number of proteins is incorrect"
assert data_prot.X.sum(axis=1).min() > 950, "Number of reads per cell is incorrect"
assert data_prot.var.Raw_Reads.sum() >= 99990 and data_prot.var.Raw_Reads.sum() <= 100010, \
  f"Expected 100000 Prot reads, got {data_prot.var.Raw_Reads.sum()}"


# check smk data
expected_sample_tags = (["SampleTag01_hs", "SampleTag02_hs", "SampleTag03_hs"] * 34)[:100]
expected_sample_names = (["Sample1", "Sample2", "Sample3"] * 34)[:100]
sample_tags = data_rna.obs["Sample_Tag"]
assert sample_tags.nunique() == 3, "Number of sample tags is incorrect"
assert sample_tags.tolist() == expected_sample_tags, "Sample tags are incorrect"
sample_names = data_rna.obs["Sample_Name"]
assert sample_names.nunique() == 3, "Number of sample names is incorrect"
assert sample_names.tolist() == expected_sample_names, "Sample names are incorrect"

# TODO: add VDJ, ATAC, and targeted RNA to test

#########################################################################################

print("> Test successful", flush=True)
