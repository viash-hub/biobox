test_dir="$meta_resources_dir/test_data"

echo "> Running busco with lineage dataset"

"$meta_executable" \
    --input $test_dir/protein.fasta \
    --mode proteins \
    --lineage_dataset stramenopiles_odb10 \
    --output_dir output 

echo ">> Checking output"
[ ! -f "output/full_table.tsv" ] && echo "full_table.tsv does not exist" && exit 1
[ ! -f "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist" && exit 1
[ ! -f "output/short_summary.json" ] && echo "short_summary.json does not exist" && exit 1
[ ! -f "output/short_summary.txt" ] && echo "short_summary.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output/full_table.tsv" ] && echo "full_table.tsv is empty" && exit 1
[ ! -s "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv is empty" && exit 1
[ ! -s "output/short_summary.json" ] && echo "short_summary.json is empty" && exit 1
[ ! -s "output/short_summary.txt" ] && echo "short_summary.txt is empty" && exit 1

echo "> Running busco with auto lineage"

"$meta_executable" \
    --input $test_dir/protein.fasta \
    --mode proteins \
    --auto_lineage \
    --output_dir output 

echo "> Finished running busco with auto lineage"
ls -l output

echo ">> Checking output"
[ ! -f "output/full_table.tsv" ] && echo "full_table.tsv does not exist" && exit 1
[ ! -f "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist" && exit 1
[ ! -f "output/short_summary.json" ] && echo "short_summary.json does not exist" && exit 1
[ ! -f "output/short_summary.txt" ] && echo "short_summary.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output/full_table.tsv" ] && echo "full_table.tsv is empty" && exit 1
[ ! -s "output/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv is empty" && exit 1
[ ! -s "output/short_summary.json" ] && echo "short_summary.json is empty" && exit 1
[ ! -s "output/short_summary.txt" ] && echo "short_summary.txt is empty" && exit 1

rm -r output/ 