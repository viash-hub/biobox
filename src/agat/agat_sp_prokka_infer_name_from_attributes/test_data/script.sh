#!/bin/bash

# Minimal prokka file
output_file="src/agat/agat_sp_prokka_infer_name_from_attributes/test_data/input_prokka.gff"

cat <<EOL > $output_file
##gff-version 3
NZ_CP006571.1	Prodigal:002006	CDS	2682	4016	.	+	0	ID=DEKMBENF_00004;eC_number=2.5.1.7;db_xref=COG:COG0766;gene=murA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P45025;locus_tag=DEKMBENF_00004;product=UDP-N-acetylglucosamine 1-carboxyvinyltransferase
EOL

# expected output: name is inferred from attributes
output_file="src/agat/agat_sp_prokka_infer_name_from_attributes/test_data/output_prokka.gff"

cat <<EOL > $output_file
##gff-version 3
NZ_CP006571.1	AGAT	gene	2682	4016	.	+	.	ID=agat-gene-1;Name=murA;db_xref=COG:COG0766;eC_number=2.5.1.7;gene=murA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P45025;locus_tag=DEKMBENF_00004;product=UDP-N-acetylglucosamine 1-carboxyvinyltransferase
NZ_CP006571.1	AGAT	mRNA	2682	4016	.	+	.	ID=agat-rna-1;Parent=agat-gene-1;db_xref=COG:COG0766;eC_number=2.5.1.7;gene=murA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P45025;locus_tag=DEKMBENF_00004;product=UDP-N-acetylglucosamine 1-carboxyvinyltransferase
NZ_CP006571.1	AGAT	exon	2682	4016	.	+	.	ID=agat-exon-1;Parent=agat-rna-1;db_xref=COG:COG0766;eC_number=2.5.1.7;gene=murA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P45025;locus_tag=DEKMBENF_00004;product=UDP-N-acetylglucosamine 1-carboxyvinyltransferase
NZ_CP006571.1	Prodigal:002006	CDS	2682	4016	.	+	0	ID=DEKMBENF_00004;Parent=agat-rna-1;db_xref=COG:COG0766;eC_number=2.5.1.7;gene=murA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P45025;locus_tag=DEKMBENF_00004;product=UDP-N-acetylglucosamine 1-carboxyvinyltransferase
EOL