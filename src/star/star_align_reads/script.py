import tempfile
import subprocess
import shutil
from pathlib import Path

## VIASH START
par = {
    "input": [
        "src/star/star_align_reads/test_data/a_R1.1.fastq",
        "src/star/star_align_reads/test_data/a_R1.2.fastq",
    ],
    "input_r2": [
        "src/star/star_align_reads/test_data/a_R2.1.fastq",
        "src/star/star_align_reads/test_data/a_R2.2.fastq",
    ],
    "genomeDir": "src/star/star_align_reads/test_data/genome.fasta",
    "aligned_reads": "aligned_reads.sam"
}
meta = {
    "cpus": 8,
    "temp_dir": "/tmp"
}
## VIASH END

##################################################
# check and process SE / PE R1 input files
input_r1 = par["input"]
readFilesIn = ",".join(par["input"])
par["input"] = None

# check and process PE R2 input files
input_r2 = par["input_r2"]
if input_r2 is not None:
    if len(input_r1) != len(input_r2):
        raise ValueError("The number of R1 and R2 files do not match.")
    readFilesIn = [readFilesIn, ",".join(par["input_r2"])]
    par["input_r2"] = None

# store readFilesIn
par["readFilesIn"] = readFilesIn

##################################################

# determine readFilesCommand
if input_r1[0].endswith(".gz"):
    print(">> Input files are gzipped, setting readFilesCommand to zcat", flush=True)
    par["readFilesCommand"] = "zcat"
elif input_r1[0].endswith(".bz2"):
    print(">> Input files are bzipped, setting readFilesCommand to bzcat", flush=True)
    par["readFilesCommand"] = "bzcat"

##################################################
# store output paths
expected_outputs = {
    "aligned_reads": ["Aligned.out.sam", "Aligned.out.bam"],
    "reads_per_gene": "ReadsPerGene.out.tab",
    "chimeric_junctions": "Chimeric.out.junction",
    "log": "Log.final.out",
    "splice_junctions": "SJ.out.tab",
    "unmapped": "Unmapped.out.mate1",
    "unmapped_r2": "Unmapped.out.mate2", 
    "reads_aligned_to_transcriptome": "Aligned.toTranscriptome.out.bam"
}
output_paths = {name: par[name] for name in expected_outputs.keys()}
for name in expected_outputs.keys():
    par[name] = None

##################################################
# process other args
par["runMode"] = "alignReads"

if "cpus" in meta and meta["cpus"]:
    par["runThreadN"] = meta["cpus"]

##################################################
# run STAR and move output to final destination
with tempfile.TemporaryDirectory(prefix="star-", dir=meta["temp_dir"], ignore_cleanup_errors=True) as temp_dir:
    print(">> Constructing command", flush=True)

    # set output paths
    temp_dir = Path(temp_dir)
    par["outTmpDir"] = temp_dir / "tempdir"
    out_dir = temp_dir / "out"
    par["outFileNamePrefix"] = f"{out_dir}/" # star needs this slash

    # construct command
    cmd_args = [ "STAR" ]
    for name, value in par.items():
        if value is not None:
            val_to_add = value if isinstance(value, list) else [value]
            cmd_args.extend([f"--{name}"] + [str(x) for x in val_to_add])
    print("", flush=True)

    # run command
    print(">> Running STAR with command:", flush=True)
    print(f"+ {' '.join(cmd_args)}", end="\n\n", flush=True)
    subprocess.run(
        cmd_args,
        check=True
    )
    print(">> STAR finished successfully", end="\n\n", flush=True)

    # move output to final destination
    print(">> Moving output to final destination", flush=True)
    for name, paths in expected_outputs.items():
        for expected_path in [paths] if isinstance(paths, str) else paths:
            expected_full_path = out_dir / expected_path
            if output_paths[name] and expected_full_path.is_file():
                print(f">> Moving {expected_path} to {output_paths[name]}", flush=True)
                shutil.move(expected_full_path, output_paths[name])
