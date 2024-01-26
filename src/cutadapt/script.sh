#!/bin/bash

if [ -z $par_outputDir ]; then
	par_outputDir=.
else
	mkdir -p "$par_outputDir"
fi

for f in $par_input; do
	[ ! -f "$f" ] && echo "The input file $f does not exist" && exit 1
done

barcodesFasta="barcodes.fasta"

awk '{print ">"$1"\n""^"$1}' $par_barcodesFile >$barcodesFasta

fastqFiles=$(echo $par_input | tr " " "\n")
for file in $fastqFiles; do
	if echo "$file" | grep -q R1; then
		input_R1=$(echo $file | grep R1)
	fi
	if echo "$file" | grep -q R2; then
		input_R2=$(echo $file | grep R2)
	fi
done
demuxFilesIn="$input_R1 $input_R2"

# Note to self:
#   The eval is here to expand shell globs, this way it is possible to use
#   for instance pointers to ".../...R?....fastq", but please use the double
#   quotes and an absolute path!
eval /usr/local/bin/cutadapt \
	-e "$par_e" \
	--no-indels \
	--action=none \
	--cores=0 \
	-g "file:$barcodesFasta" \
	-o "$par_outputDir/{name}_R1_001.fastq" \
	-p "$par_outputDir/{name}_R2_001.fastq" \
	"$demuxFilesIn" >"$par_report"
