set -e

test_dir="$meta_resources_dir/test_data"

mkdir "run_qualimap_rnaseq_html"
cd "run_qualimap_rnaseq_html"

echo "> Running qualimap with html output report"

"$meta_executable" \
    --input $test_dir/a.bam \
    --gtf $test_dir/annotation.gtf \
    --output_report report.html \
    --output_counts counts.txt \
    --output output.txt

echo ">> Checking output"
[ ! -f "report.html" ] && echo "report.html does not exist" && exit 1
[ ! -f "counts.txt" ] && echo "counts.txt does not exist" && exit 1
[ ! -f "output.txt" ] && echo "output.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "report.html" ] && echo "report.html is empty" && exit 1
[ ! -s "counts.txt" ] && echo "counts.txt is empty" && exit 1
[ ! -s "output.txt" ] && echo "output.txt is empty" && exit 1

cd ..
rm -r run_qualimap_rnaseq_html

mkdir "run_qualimap_rnaseq_pdf"
cd "run_qualimap_rnaseq_pdf"

echo "> Running qualimap with pdf output report"

"$meta_executable" \
    --input $test_dir/a.bam \
    --gtf $test_dir/annotation.gtf \
    --output_report report.pdf \
    --output_counts counts.txt \
    --output output.txt

echo ">> Checking output"
[ ! -f "report.pdf" ] && echo "report.pdf does not exist" && exit 1
[ ! -f "counts.txt" ] && echo "counts.txt does not exist" && exit 1
[ ! -f "output.txt" ] && echo "output.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "report.pdf" ] && echo "report.pdf is empty" && exit 1
[ ! -s "counts.txt" ] && echo "counts.txt is empty" && exit 1
[ ! -s "output.txt" ] && echo "output.txt is empty" && exit 1

cd ..
rm -r run_qualimap_rnaseq_pdf

mkdir "run_qualimap_rnaseq"
cd "run_qualimap_rnaseq"

echo "> Running qualimap without report and counts output"

"$meta_executable" \
    --input $test_dir/a.bam \
    --gtf $test_dir/annotation.gtf \
    --output output.txt

echo ">> Checking output"
[ -f "report.pdf" ] && echo "report.pdf exists" && exit 1
[ -f "counts.txt" ] && echo "counts.txt exists" && exit 1
[ ! -f "output.txt" ] && echo "output.txt does not exist" && exit 1

echo ">> Checking if output is empty"
[ ! -s "output.txt" ] && echo "output.txt is empty" && exit 1

cd ..
rm -r run_qualimap_rnaseq