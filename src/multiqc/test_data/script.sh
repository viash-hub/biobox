# multiqc test data

# Test data was obtained from https://multiqc.info/example-reports/

curl -O -J -L http://multiqc.info/examples/rna-seq/data.zip
unzip -q "data.zip" -d "src/multiqc/test_data/rna-seq"

