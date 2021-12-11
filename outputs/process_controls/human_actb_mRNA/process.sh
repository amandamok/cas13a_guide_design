#!/bin/bash

rm *txt
rm *fa
rm *sam
rm *bowtiestats

scripts_dir="../../../scripts"
ref_dir="../../../ref_data"

# generate spacer windows
Rscript $scripts_dir/break_genome.R -i mature_mRNA.fasta -w 20 -e Cas13a -g + -s +

# score specificity against human coronaviruses
Rscript $scripts_dir/score_human_CoV_specificity.R -e Cas13a

# align against human transcriptome
grep "(ACTB)" $ref_dir/GRCh38_latest_rna.fna | cut -f1 -d " " | cut -f2 -d ">" > GRCh38_actb.txt
Rscript $scripts_dir/align_bowtie.R -g GRCh38_latest_rna -m 2 -e Cas13a -v GRCh38_actb.txt

# align against cow transcriptome
Rscript $scripts_dir/align_bowtie.R -g ARS-UCD1_rna -m 2 -e Cas13a

# calculate crRNA folding structures (RNAfold)
Rscript $scripts_dir/score_RNAfold_crRNAs.R -e Cas13a

# compile scores and add extra annotations

