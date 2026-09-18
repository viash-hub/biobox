import csv
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

## VIASH START
par = {
    "csv": "config.csv",
    "fastqs": ["cellranger_tiny_fastq"],
    "gex_reference": "cellranger_tiny_ref",
    "probe_set": None,
    "cmo_set": None,
    "barcode_sample_assignment": None,
    "tenx_cloud_token": None,
    "feature_reference": None,
    "vdj_reference": None,
    "inner_enrichment_primers": None,
    "output": "output",
    "id": "run",
    "description": None,
    "dry": False,
}
meta = {
    "cpus": 4,
    "memory_gb": 8,
    "temp_dir": "/tmp",
    "name": "cellranger_multi",
}
## VIASH END

# The multi config CSV refers to its inputs by path. Every entry listed here is rewritten to
# the staged location of the corresponding component argument. The alias is the alternative
# spelling that Cell Ranger accepts for the same entry.
STAGED_ENTRIES = [
    # (argument, section, entry, alias)
    ("gex_reference", "gene-expression", "reference", "ref"),
    ("probe_set", "gene-expression", "probe-set", None),
    ("cmo_set", "gene-expression", "cmo-set", None),
    ("barcode_sample_assignment", "gene-expression", "barcode-sample-assignment", None),
    ("tenx_cloud_token", "gene-expression", "tenx-cloud-token-path", None),
    ("feature_reference", "feature", "reference", "ref"),
    ("vdj_reference", "vdj", "reference", "ref"),
    ("inner_enrichment_primers", "vdj", "inner-enrichment-primers", None),
]

# These are passed to Cell Ranger as a directory, so they cannot be used while packed.
UNPACKED_ARGUMENTS = ["gex_reference", "vdj_reference"]

# The file names Cell Ranger accepts are `<fastq_id>_S<n>[_L<lane>]_<R1|R2|R3|RA|I1|I2>_
# <chunk>.fastq`, optionally `.gz` or `.lz4` compressed. The `.fq` extension is not accepted.
FASTQ_PATTERNS = ["*.fastq", "*.fastq.gz", "*.fastq.lz4"]

SECTION_HEADER = re.compile(r"^\[(.+)\]$")

# Entries of the multi config CSV, keyed by `(section, entry)`.
SectionMap = dict[tuple[str, str], str]


def resolve_reference(reference: str, name: str, workdir: Path) -> Path:
    """Resolve a reference to the directory Cell Ranger expects.

    References are directories, but are commonly distributed as a `.tar.gz` archive. Those
    are unpacked into the working directory. Such an archive normally holds a single top
    level directory, which is the reference itself.
    """
    path = Path(reference).resolve()
    if path.is_dir() or not tarfile.is_tarfile(path):
        return path

    print(f"> Untarring {name}", flush=True)
    unpack_dir = workdir / name
    unpack_dir.mkdir(parents=True)
    with tarfile.open(path) as archive:
        archive.extractall(unpack_dir, filter="data")

    entries = list(unpack_dir.iterdir())
    if len(entries) == 1 and entries[0].is_dir():
        return entries[0]
    return unpack_dir


def stage_fastqs(fastqs: list[str], fastq_dir: Path) -> dict[Path, Path]:
    """Symlink the FASTQ files into one staged directory per directory they come from.

    Cell Ranger identifies a set of reads by the directory that holds it and by the file
    names in it, and the same library sequenced on more than one flowcell yields identically
    named files. Only the `fastqs` column of the `[libraries]` section tells those apart, so
    the grouping of the input has to survive staging. The result maps each directory the
    files came from onto the staged directory that holds them, which `resolve_fastq_group`
    matches the column against. Directories are told apart by their full path, so the ones
    that share a name, like `flowcell1/a` and `flowcell2/a`, stay separate.
    """
    groups: dict[Path, Path] = {}
    for fastq in fastqs:
        path = Path(fastq).resolve()
        if path.is_dir():
            found = sorted(match for p in FASTQ_PATTERNS for match in path.rglob(p))
        else:
            found = [path]
        if not found:
            raise ValueError(f"No FASTQ files were found in '{fastq}'.")
        for source in found:
            group = groups.get(source.parent)
            if group is None:
                # The staged name only has to be unique; the column is matched against the
                # original path. Keep it recognisable in the config and the log all the same.
                taken = {staged.name for staged in groups.values()}
                name = source.parent.name
                while name in taken:
                    name = f"{source.parent.name}_{len(taken) + 1}"
                    taken.add(source.parent.name)
                group = groups[source.parent] = fastq_dir / name
                group.mkdir(parents=True)
            link = group / source.name
            # The same file can be passed twice, as a directory and as a file in it.
            if not link.exists():
                link.symlink_to(source)
    return groups


def shared_tail(parts: tuple[str, ...], other: tuple[str, ...]) -> int:
    """Count the trailing path components that two paths have in common."""
    shared = 0
    while (
        shared < len(parts)
        and shared < len(other)
        and parts[-1 - shared] == other[-1 - shared]
    ):
        shared += 1
    return shared


def resolve_fastq_group(value: str, groups: dict[Path, Path]) -> Path:
    """Map a value of the `fastqs` column onto the directory its reads were staged in.

    The value is matched against the directories the files came from by their longest common
    trailing path, so the column only has to carry enough of the path to tell those apart,
    rather than the path on the machine that runs this. A CSV that has a placeholder there
    works as long as a single directory of FASTQ files was passed to the component.
    """
    wanted = Path(value.strip()).parts
    shared = {source: shared_tail(source.parts, wanted) for source in groups}
    longest = max(shared.values(), default=0)
    matched = [source for source, common in shared.items() if common == longest]

    if longest > 0 and len(matched) == 1:
        return groups[matched[0]]
    if longest == 0 and len(groups) == 1:
        return next(iter(groups.values()))
    raise ValueError(
        f"The fastqs column of the [libraries] section refers to '{value}', which does not "
        "match exactly one of the directories the FASTQ files came from "
        f"({', '.join(sorted(str(source) for source in groups))}). Point it at the directory "
        "that holds the reads of that library, with enough of the path to tell it apart."
    )


def parse_sections(rows) -> tuple[list[list[str]], list[tuple[str, list[list[str]]]]]:
    """Split the rows of a multi config CSV into the rows before the first section header and
    the sections themselves, as `(name, rows)` pairs.

    Blank lines are not significant to Cell Ranger, so they are dropped here and one is
    written back out between the sections. Section names are not case sensitive either.
    """
    preamble: list[list[str]] = []
    sections: list[tuple[str, list[list[str]]]] = []
    for row in rows:
        if not any(field.strip() for field in row):
            continue
        header = SECTION_HEADER.match(row[0].strip())
        if header:
            sections.append((header.group(1).lower(), []))
        elif sections:
            sections[-1][1].append(row)
        else:
            preamble.append(row)
    return preamble, sections


def take_entries(section: str, pending: SectionMap) -> list[list[str]]:
    """Remove the entries that are still pending for a section and format them as rows."""
    return [
        [entry, pending.pop((section, entry))]
        for name, entry in list(pending)
        if name == section
    ]


def stage_entries(
    section: str, rows: list[list[str]], pending: SectionMap, aliases: SectionMap
) -> None:
    """Point the entries of an option section at the staged input files, in place.

    An entry that the section does not have is appended to it. A comment never matches an
    entry that was staged, so it passes through like any other row.
    """
    for index, row in enumerate(rows):
        entry = row[0].strip().lower()
        entry = aliases.get((section, entry), entry)
        if (section, entry) in pending:
            rows[index] = [row[0].strip(), pending.pop((section, entry))]
    rows += take_entries(section, pending)


def stage_libraries(rows: list[list[str]], fastq_groups: dict[Path, Path]) -> None:
    """Point every value of the `fastqs` column at the directory its reads were staged in."""
    fastq_column = None
    for row in rows:
        if row[0].lstrip().startswith("#"):
            continue
        if fastq_column is None:
            # The first row of the section holds the column names.
            columns = [field.strip().lower() for field in row]
            if "fastqs" not in columns:
                raise ValueError(
                    "The [libraries] section of the multi config CSV has no fastqs column, "
                    f"found: {', '.join(columns)}."
                )
            fastq_column = columns.index("fastqs")
            continue
        row[fastq_column] = str(resolve_fastq_group(row[fastq_column], fastq_groups))


def rewrite_multi_config(
    config: str,
    overrides: list[tuple[str, str, str | None, str]],
    fastq_groups: dict[Path, Path],
    destination: Path,
) -> None:
    """Write a copy of the multi config CSV that points at the staged input files.

    Every override, a `(section, entry, alias, value)` tuple, replaces the value of that
    entry, and every value of the `fastqs` column of the `[libraries]` section is replaced by
    the staged directory it names. An entry that its section does not have is appended to
    that section, and a section that the CSV does not have is appended to the CSV.
    """
    with open(config, newline="", encoding="utf-8") as handle:
        preamble, sections = parse_sections(csv.reader(handle))

    pending = {(section, entry): value for section, entry, _, value in overrides}
    aliases = {
        (section, alias): entry for section, entry, alias, _ in overrides if alias is not None
    }

    for name, rows in sections:
        if name == "libraries":
            stage_libraries(rows, fastq_groups)
        else:
            stage_entries(name, rows, pending, aliases)

    # Whatever is left belongs to a section that the CSV does not have at all.
    for name in dict.fromkeys(section for section, _ in pending):
        sections.append((name, take_entries(name, pending)))

    with open(destination, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(preamble)
        for name, rows in sections:
            writer.writerow([f"[{name}]"])
            writer.writerows(rows)
            writer.writerow([])


with tempfile.TemporaryDirectory(dir=meta["temp_dir"], prefix=f"{meta['name']}-") as temp_dir:
    workdir = Path(temp_dir)

    for argument in UNPACKED_ARGUMENTS:
        if par[argument] is not None:
            par[argument] = resolve_reference(par[argument], argument, workdir)

    overrides = [
        (section, entry, alias, str(Path(par[argument]).resolve()))
        for argument, section, entry, alias in STAGED_ENTRIES
        if par[argument] is not None
    ]

    fastq_groups = stage_fastqs(par["fastqs"], workdir / "fastqs")

    print("> Staging input files in the multi config CSV", flush=True)
    config_csv = workdir / "config.csv"
    rewrite_multi_config(par["csv"], overrides, fastq_groups, config_csv)

    cmd = [
        "cellranger",
        "multi",
        "--id",
        par["id"],
        "--csv",
        str(config_csv),
        "--disable-ui",
    ]
    if meta["cpus"]:
        cmd += ["--localcores", str(meta["cpus"])]
    # Cell Ranger needs some headroom on top of the memory it is allowed to request.
    if meta["memory_gb"]:
        if meta["memory_gb"] < 2:
            print("WARNING: Memory is less than 2GB, unsetting memory requirements", flush=True)
        else:
            cmd += ["--localmem", str(meta["memory_gb"] - 2)]
    if par["description"] is not None:
        cmd += ["--description", par["description"]]
    if par["dry"]:
        cmd.append("--dry")

    print("> Running cellranger multi", flush=True)
    output = Path(par["output"])
    output.mkdir(parents=True, exist_ok=True)
    # Cell Ranger keeps its own log in the pipestance directory, which is inside the
    # temporary directory and is discarded with it, so stream a copy next to the results.
    # The log is written while the pipeline runs, so it survives a failed run as well.
    with (
        (output / "cellranger_multi.log").open("w", buffering=1) as open_log,
        subprocess.Popen(
            cmd,
            cwd=workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            errors="replace",  # Handle encoding errors in stderr and stdout
            encoding="utf-8",
        ) as process,
    ):
        for line in process.stdout:
            print(line, end="", flush=True)
            open_log.write(line)
    if process.returncode != 0:
        raise RuntimeError(
            f"cellranger multi returned a nonzero exitcode ({process.returncode})."
        )

    print("> Copying output", flush=True)
    outs = workdir / par["id"] / "outs"
    if par["dry"]:
        # The pipeline was not run, so there is no output to copy. Hand over the multi config
        # CSV that was generated for it, which is what Cell Ranger copies there otherwise.
        shutil.copy(config_csv, output / "config.csv")
    else:
        if not outs.is_dir():
            raise RuntimeError(
                f"Cell Ranger reported success but produced no output directory at '{outs}'."
            )
        for item in outs.iterdir():
            shutil.move(item, output)
