set -e

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_doesnt_exist() {
  [ ! -f "$1" ] || { echo "File '$1' exists but shouldn't" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
#############################################


test_dir="$meta_resources_dir/test_data"

mkdir "run_qualimap_rnaseq_html"
cd "run_qualimap_rnaseq_html"

echo "> Running qualimap with html output report"

"$meta_executable" \
    --bam $test_dir/a.bam \
    --gtf $test_dir/annotation.gtf \
    --report report.html \
    --counts counts.txt \
    --qc_results output.txt

echo ">> Checking output"
assert_file_exists "report.html"
assert_file_exists "counts.txt"
assert_file_exists "output.txt"
assert_file_doesnt_exist "report.pdf"

echo ">> Checking if output is empty"
assert_file_not_empty "report.html"
assert_file_not_empty "counts.txt"
assert_file_not_empty "output.txt"

echo ">> Checking output contents"
assert_file_contains "output.txt" ">>>>>>> Input"
assert_file_contains "output.txt" ">>>>>>> Reads alignment"
assert_file_contains "output.txt" ">>>>>>> Reads genomic origin"
assert_file_contains "output.txt" ">>>>>>> Transcript coverage profile"
assert_file_contains "output.txt" ">>>>>>> Junction analysis"
assert_file_contains "output.txt" ">>>>>>> Transcript coverage profile"

assert_file_contains "counts.txt" "ENSG00000125841.12"

assert_file_contains "report.html" "<title>Qualimap report: RNA Seq QC</title>"
assert_file_contains "report.html" "<h3>Input</h3>"
assert_file_contains "report.html" "<h3>Reads alignment</h3>"
assert_file_contains "report.html" "<h3>Reads genomic origin</h3>"
assert_file_contains "report.html" "<h3>Transcript coverage profile</h3>"
assert_file_contains "report.html" "<h3>Junction analysis</h3>"


cd ..
rm -r run_qualimap_rnaseq_html

mkdir "run_qualimap_rnaseq_pdf"
cd "run_qualimap_rnaseq_pdf"

echo "> Running qualimap with pdf output report"

"$meta_executable" \
    --bam $test_dir/a.bam \
    --gtf $test_dir/annotation.gtf \
    --report report.pdf \
    --counts counts.txt \
    --qc_results output.txt

echo ">> Checking output"
assert_file_exists "report.pdf"
assert_file_exists "counts.txt"
assert_file_exists "output.txt"
assert_file_doesnt_exist "report.html"

echo ">> Checking if output is empty"
assert_file_not_empty "report.pdf"
assert_file_not_empty "counts.txt"
assert_file_not_empty "output.txt"

cd ..
rm -r run_qualimap_rnaseq_pdf

mkdir "run_qualimap_rnaseq"
cd "run_qualimap_rnaseq"

echo "> Running qualimap without report and counts output"

"$meta_executable" \
    --bam $test_dir/a.bam \
    --gtf $test_dir/annotation.gtf \
    --qc_results output.txt

echo ">> Checking output"
assert_file_doesnt_exist "report.pdf"
assert_file_doesnt_exist "report.html"
assert_file_doesnt_exist "counts.txt"
assert_file_exists "output.txt"

echo ">> Checking if output is empty"
assert_file_not_empty "output.txt"

cd ..
rm -r run_qualimap_rnaseq