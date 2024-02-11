echo "> Downloading busco stramenopiles_odb10 dataset"

"$meta_executable" \
    --download stramenopiles_odb10 \
    --download_path downloads 

echo ">> Checking output"
[ ! -f "downloads/file_versions.tsv" ] && echo "file_versions.tsv does not exist" && exit 1
[ ! -f "downloads/lineages/stramenopiles_odb10/dataset.cfg" ] && echo "dataset.cfg does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "downloads/file_versions.tsv" ] && echo "file_versions.tsv is empty" && exit 1
[ ! -s "downloads/lineages/stramenopiles_odb10/dataset.cfg" ] && echo "dataset.cfg is empty" && exit 1

rm -r downloads