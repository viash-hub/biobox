test_dir="$meta_resources_dir/test_data"

echo "> Running busco"
"$meta_executable" \
    --input $test_dir/protein.fasta \
    --mode protein \
    --lineage_dataset stramenopiles_odb10 \
    --output_dir output 

echo ">> Checking output"
[ ! -f "output/short_summary.specific.stramenopiles_odb10.protein.fasta.json" ] && echo "specific_short_summary.json does not exist" && exit 1
[ ! -f "output/short_summary.specific.stramenopiles_odb10.protein.fasta.txt" ] && echo "specific_short_summary.txt does not exist" && exit 1
[ ! -f "output/run_stramenopiles_odb10/full_table.tsv" ] && echo "full_table.tsv does not exist" && exit 1
[ ! -f "output/run_stramenopiles_odb10/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv does not exist" && exit 1
[ ! -f "output/run_stramenopiles_odb10/short_summary.json" ] && echo "short_summary.json does not exist" && exit 1
[ ! -f "output/run_stramenopiles_odb10/short_summary.txt" ] && echo "short_summary.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output/short_summary.specific.stramenopiles_odb10.protein.fasta.json" ] && echo "specific_short_summary.json is empty" && exit 1
[ ! -s "output/short_summary.specific.stramenopiles_odb10.protein.fasta.txt" ] && echo "specific_short_summary.txt is empty" && exit 1
[ ! -s "output/run_stramenopiles_odb10/full_table.tsv" ] && echo "full_table.tsv is empty" && exit 1
[ ! -s "output/run_stramenopiles_odb10/missing_busco_list.tsv" ] && echo "missing_busco_list.tsv is empty" && exit 1
[ ! -s "output/run_stramenopiles_odb10/short_summary.json" ] && echo "short_summary.json is empty" && exit 1
[ ! -s "output/run_stramenopiles_odb10/short_summary.txt" ] && echo "short_summary.txt is empty" && exit 1


rm -r output/ 
