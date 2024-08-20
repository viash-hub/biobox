import os
import re
import subprocess
import tempfile
from typing import Any
import yaml
import shutil
import glob

## VIASH START
par = {
    'reads': [
        'resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R1_001_subset.fastq.gz', 
        'resources_test/bdrhap_5kjrt/raw/12ABC_S1_L432_R2_001_subset.fastq.gz'
    ],
    'reads_atac': None,
    'reference_archive': "resources_test/reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz",
    'targeted_reference': [],
    'abseq_reference': [],
    'supplemental_reference': [],
    'output': 'output_dir',
    'cell_calling_data': None,
    'cell_calling_bioproduct_algorithm': None,
    'cell_calling_atac_algorithm': None,
    'exact_cell_count': None,
    'expected_cell_count': None,
    'exclude_intronic_reads': None,
    'sample_tags_version': None,
    'tag_names': [],
    'vdj_version': None,
    'predefined_atac_peaks': None,
    'run_name': "sample",
    'generate_bam': None,
    'alignment_star_params': None,
    'alignment_bwa_mem2_params': None,
    'parallel': True,
    'timestamps': False,
    'dryrun': False
}
meta = {
    'config': "target/nextflow/bd_rhaspody/bd_rhaspody_sequence_analysis/.config.vsh.yaml",
    'resources_dir': os.path.abspath('src/bd_rhaspody/bd_rhaspody_sequence_analysis'),
    'temp_dir': os.getenv("VIASH_TEMP"),
    'memory_mb': None,
    'cpus': None
}
## VIASH END

def clean_arg(argument):
    argument["clean_name"] = argument["name"].lstrip("-")
    return argument

def read_config(path: str) -> dict[str, Any]:
    with open(path, 'r') as f:
        config = yaml.safe_load(f)
    
    config["arguments"] = [
        clean_arg(arg)
        for grp in config["argument_groups"]
        for arg in grp["arguments"]
    ]
    
    return config

def strip_margin(text: str) -> str:
    return re.sub('(\n?)[ \t]*\|', '\\1', text)

def process_params(par: dict[str, Any], config, temp_dir: str) -> str:
    # check input parameters
    assert par["reads"] or par["reads_atac"], "Pass at least one set of inputs to --reads or --reads_atac."

    # output to temp dir if output_dir was not passed
    if not par["output_dir"]:
        par["output_dir"] = os.path.join(temp_dir, "output")

    # checking sample prefix
    if par["run_name"] and re.match("[^A-Za-z0-9]", par["run_name"]):
        print("--run_name should only consist of letters, numbers or hyphens. Replacing all '[^A-Za-z0-9]' with '-'.", flush=True)
        par["run_name"] = re.sub("[^A-Za-z0-9\\-]", "-", par["run_name"])

    # make paths absolute
   for argument in config["arguments"]:
        arg_clean_name = argument["clean_name"]
        if not par[arg_clean_name] or not argument["type"] == "file":
            continue
        par_value = par[arg_clean_name]
        if isinstance(par_value, list):
            par_value_absolute = list(map(os.path.abspath, file_path))
        else:
            par_value_absolute = os.path.abspath(file_path)
        par_value_absolute[arg_clean_name] = par_value_absolute
    
    return par

def generate_config(par: dict[str, Any], config) -> str:
    content_list = [strip_margin(f"""\
        |#!/usr/bin/env cwl-runner
        |
        |cwl:tool: rhapsody
        |""")]

    for argument in config["arguments"]:
        arg_clean_name = argument["clean_name"]
        arg_par_value = par[arg_clean_name]
        arg_info = argument.get("info") or {}
        config_key = arg_info.get("config_key")
        if arg_par_value and config_key:

            if argument["type"] == "file":
                str = strip_margin(f"""\
                    |{config_key}:
                    |""")
                if isinstance(arg_par_value, list):
                    for file in arg_par_value:
                        str += strip_margin(f"""\
                            | - class: File
                            |   location: "{file}"
                            |""")
                else:
                    str += strip_margin(f"""\
                        |   class: File
                        |   location: "{arg_par_value}"
                        |""")
                content_list.append(str)
            else:
                content_list.append(strip_margin(f"""\
                    |{config_key}: {arg_par_value}
                    |"""))

    ## Write config to file
    return ''.join(content_list)

def generate_config_file(par: dict[str, Any], config: dict[str, Any], temp_dir: str) -> str:
    config_file = os.path.join(temp_dir, "config.yml")
    config_content = generate_config(par, config)
    with open(config_file, "w") as f:
        f.write(config_content)
    return config_file

def generate_cwl_file(meta: dict[str, Any], dir: str) -> str:
    # create cwl file (if need be)
    orig_cwl_file=os.path.join(meta["resources_dir"], "rhapsody_pipeline_2.2.1_nodocker.cwl")
    if not meta["memory_mb"] and not meta["cpus"]:
        return os.path.abspath(orig_cwl_file)
     
    # Inject computational requirements into pipeline
    cwl_file = os.path.join(dir, "pipeline.cwl")

    # Read in the file
    with open(orig_cwl_file, 'r') as file :
        cwl_data = file.read()

    # Inject computational requirements into pipeline
    if meta["memory_mb"]:
        memory = int(meta["memory_mb"]) - 2000 # keep 2gb for OS
        cwl_data = re.sub('"ramMin": [^\n]*[^,](,?)\n', f'"ramMin": {memory}\\1\n', cwl_data)
    if meta["cpus"]:
        cwl_data = re.sub('"coresMin": [^\n]*[^,](,?)\n', f'"coresMin": {meta["cpus"]}\\1\n', cwl_data)

    # Write the file out again
    with open(cwl_file, 'w') as file:
        file.write(cwl_data)
        
    return os.path.abspath(cwl_file)

def copy_outputs(par: dict[str, Any], config: dict[str, Any]):
    for arg in config["arguments"]:
        par_value = par[arg["clean_name"]]
        if par_value and arg["type"] == "file" and arg["direction"] == "output":
            # example template: '[sample_name]_(assay)_cell_type_experimental.csv'
            template = (arg.get("info") or {}).get("template")
            if template:
                template_glob = template\
                    .replace("[sample_name]", par["run_name"])\
                    .replace("(assay)", "*")\
                    .replace("[number]", "*")
                files = glob.glob(os.path.join(par["output_dir"], template_glob))
                if not files and arg["required"]:
                    raise ValueError(f"Expected output file '{template_glob}' not found.")
                elif len(files) > 1 and not arg["multiple"]:
                    raise ValueError(f"Expected single output file '{template_glob}', but found multiple.")
                
                if not arg["multiple"]:
                    shutil.copy(files[0], par_value)
                else:
                    # replace '*' in par_value with index
                    for i, file in enumerate(files):
                        shutil.copy(file, par_value.replace("*", str(i)))


def main(par: dict[str, Any], meta: dict[str, Any], temp_dir: str):
    config = read_config(meta["config"])
    
    # Preprocess params
    par = process_params(par, config, temp_dir)

    ## Process parameters
    cmd = [
        "cwl-runner",
        "--no-container",
        "--preserve-entire-environment",
        "--outdir", par["output_dir"],
    ]

    if par["parallel"]:
        cmd.append("--parallel")

    if par["timestamps"]:
        cmd.append("--timestamps")

    # Create cwl file (if need be)
    cwl_file = generate_cwl_file(meta, temp_dir)
    cmd.append(cwl_file)

    # Create params file
    config_file = generate_config_file(par, config, temp_dir)
    cmd.append(config_file)
    
    # keep environment variables but set TMPDIR to temp_dir
    env = dict(os.environ)
    env["TMPDIR"] = temp_dir

    # Create output dir if not exists
    if not os.path.exists(par["output_dir"]):
        os.makedirs(par["output_dir"])

    # Run command
    print("> " + ' '.join(cmd), flush=True)
    _ = subprocess.check_call(
        cmd,
        cwd=os.path.dirname(config_file),
        env=env
    )

    # Copy outputs
    copy_outputs(par, config)

if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="cwl-bd_rhapsody-", dir=meta["temp_dir"]) as temp_dir:
        main(par, meta, temp_dir)
