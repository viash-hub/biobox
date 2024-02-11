import tempfile
import subprocess
import shutil
from pathlib import Path

## VIASH START
par = {
    'input': [
        'src/star/star_align/test_data/a_R1.1.fastq',
        'src/star/star_align/test_data/a_R1.2.fastq',
    ],
    'input_r2': [
        'src/star/star_align/test_data/a_R2.1.fastq',
        'src/star/star_align/test_data/a_R2.2.fastq',
    ],
    'genome_dir': 'src/star/star_align/test_data/genome.fasta',
    'output': 'test_output'
}
meta = {
    'cpus': 8,
    'temp_dir': '/tmp'
}
## VIASH END

# check and process SE / PE R1 input files
fq1 = par["input"]
input_str = ','.join(par["input"])
par["input"] = None

# check and process PE R2 input files
fq2 = par["input_r2"]
if fq2 is not None:
    if len(fq1) != len(fq2):
        raise ValueError("The number of R1 and R2 files do not match.")
    input_str = [input_str, ','.join(par["input_r2"])]
    par["input_r2"] = None

# store readFilesIn
par["readFilesIn"] = input_str

# determine readFilesCommand
if fq1[0].endswith(".gz"):
    print(">> Input files are gzipped, setting readFilesCommand to zcat", flush=True)
    par["readFilesCommand"] = "zcat"
elif fq1[0].endswith(".bz2"):
    print(">> Input files are bzipped, setting readFilesCommand to bzcat", flush=True)
    par["readFilesCommand"] = "bzcat"

# store output paths
expected_outputs = {
    "aligned_reads": ["Aligned.out.sam", "Aligned.out.bam"],
    "reads_per_gene": "ReadsPerGene.out.tab",
    "chimeric_junctions": "Chimeric.out.junction",
    "log": "Log.final.out",
    "splice_junctions": "SJ.out.tab",
    "unmapped": "Unmapped.out.mate1",
    "unmapped_r2": "Unmapped.out.mate2"
}
output_paths = {name: par[name] for name in expected_outputs.keys()}
for name in expected_outputs.keys():
    par[name] = None

with tempfile.TemporaryDirectory(prefix="star-", dir=meta["temp_dir"], ignore_cleanup_errors=True) as temp_dir:
    print(">> Constructing command", flush=True)

    par["runMode"] = "alignReads"
    par["outTmpDir"] = temp_dir + "/tempdir"
    par["outFileNamePrefix"] = temp_dir + "/out/"
    out_dir = Path(temp_dir) / "out"

    if 'cpus' in meta and meta['cpus']:
        par["runThreadN"] = meta["cpus"]

    cmd_args = [ "STAR" ]
    for name, value in par.items():
        if value is not None:
            if isinstance(value, list):
                cmd_args.extend(["--" + name] + [str(x) for x in value])
            else:
                cmd_args.extend(["--" + name, str(value)])
    print("", flush=True)

    print(">> Running STAR with command:", flush=True)
    print("+ " + ' '.join([str(x) for x in cmd_args]), flush=True)
    print("", flush=True)

    subprocess.run(
        cmd_args,
        check=True
    )

    print(">> STAR finished successfully", flush=True)
    print("", flush=True)

    print(">> Moving output to final destination", flush=True)
    for name, paths in expected_outputs.items():
        for expected_path in [paths] if isinstance(paths, str) else paths:
            expected_full_path = out_dir / expected_path
            if output_paths[name] and expected_full_path.is_file():
                print(">> Moving " + expected_path + " to " + output_paths[name], flush=True)
                shutil.move(expected_full_path, output_paths[name])