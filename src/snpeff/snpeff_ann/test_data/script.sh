# Test files from SnpEff examples
if [ ! -f snpEff_latest_core.zip ]; then
    wget https://snpeff.odsp.astrazeneca.com/versions/snpEff_latest_core.zip
fi

if [ ! -d snpEff ]; then
    unzip snpEff_latest_core.zip
fi

mv snpEff/examples/test.vcf src/snpeff/test_data/
mv snpEff/examples/cancer.vcf src/snpeff/test_data/
mv snpEff/examples/my_annotations.bed src/snpeff/test_data/

rm -rf snpEff_latest_core.zip
rm -rf snpEff