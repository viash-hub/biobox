test_dir="$meta_resources_dir/test_data"

mkdir "run_prot_stramenopiles"
cd "run_prot_stramenopiles"

echo "> Running busco with lineage dataset"

"$meta_executable" \
    --input $test_dir/protein.fasta \
    --mode proteins \
    --lineage_dataset stramenopiles_odb10 \
    --output_dir output \
    --short_summary_json short_summary.json \
    --short_summary_txt short_summary.txt \
    --full_table full_table.tsv \
    --missing_busco_list missing_busco_list.tsv

echo ">> Checking output"
[ ! -f "output/full_table.tsv" ] && echo "full_table.tsv does not exist" && exit 1
[ ! -f "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist" && exit 1
[ ! -f "output/short_summary.json" ] && echo "short_summary.json does not exist" && exit 1
[ ! -f "output/short_summary.txt" ] && echo "short_summary.txt does not exist" && exit 1
[ ! -f "full_table.tsv" ] && echo "full_table.tsv does not exist" && exit 1
[ ! -f "missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist" && exit 1
[ ! -f "short_summary.json" ] && echo "short_summary.json does not exist" && exit 1
[ ! -f "short_summary.txt" ] && echo "short_summary.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output/full_table.tsv" ] && echo "full_table.tsv is empty" && exit 1
[ ! -s "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv is empty" && exit 1
[ ! -s "output/short_summary.json" ] && echo "short_summary.json is empty" && exit 1
[ ! -s "output/short_summary.txt" ] && echo "short_summary.txt is empty" && exit 1
[ ! -s "full_table.tsv" ] && echo "full_table.tsv is empty" && exit 1
[ ! -s "missing_busco_list.tsv" ] && echo "missing_busco_list.tsv is empty" && exit 1
[ ! -s "short_summary.json" ] && echo "short_summary.json is empty" && exit 1
[ ! -s "short_summary.txt" ] && echo "short_summary.txt is empty" && exit 1

cd ..
mkdir "run_prot_autolineage"
cd "run_prot_autolineage"

echo "> Running busco with auto lineage"

"$meta_executable" \
    --input $test_dir/protein.fasta \
    --mode proteins \
    --auto_lineage \
    --output_dir output 

echo ">> Checking output"
[ ! -f "output/full_table.tsv" ] && echo "full_table.tsv does not exist in output folder" && exit 1
[ ! -f "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist in output folder" && exit 1
[ ! -f "output/short_summary.json" ] && echo "short_summary.json does not exist in output folder" && exit 1
[ ! -f "output/short_summary.txt" ] && echo "short_summary.txt does not exist in output folder" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output/full_table.tsv" ] && echo "full_table.tsv in output folder is empty" && exit 1
[ ! -s "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv in output folder is empty" && exit 1
[ ! -s "output/short_summary.json" ] && echo "short_summary.json in output folder is empty" && exit 1
[ ! -s "output/short_summary.txt" ] && echo "short_summary.txt in output folder is empty" && exit 1

rm -r output/ 

cd ..
mkdir "run_genome"
cd "run_genome"

echo "> Running busco with genome data"

"$meta_executable" \
    --input $test_dir/genome.fna \
    --mode genome \
    --lineage_dataset saccharomycetes_odb10 \
    --output_dir output 

echo ">> Checking output"
[ ! -f "output/full_table.tsv" ] && echo "full_table.tsv does not exist in output folder" && exit 1
[ ! -f "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist in output folder" && exit 1
[ ! -f "output/short_summary.json" ] && echo "short_summary.json does not exist in output folder" && exit 1
[ ! -f "output/short_summary.txt" ] && echo "short_summary.txt does not exist in output folder" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output/full_table.tsv" ] && echo "full_table.tsv in output folder is empty" && exit 1
[ ! -s "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv in output folder is empty" && exit 1
[ ! -s "output/short_summary.json" ] && echo "short_summary.json in output folder is empty" && exit 1
[ ! -s "output/short_summary.txt" ] && echo "short_summary.txt in output folder is empty" && exit 1

rm -r output/ 