#!/bin/bash

rm *txt
rm *fa
rm *sam
rm *bowtiestats

# generate spacer windows
Rscript ../../scripts/break_genome.R -i human_pop1.fasta -w 20 -e Cas13a -g + -s + 

# score specificity against human coronaviruses
Rscript ../../scripts/score_human_CoV_specificity.R

# align against human transcriptome
grep POP1 ../../ref_data/GRCh38_latest_rna.fna | cut -f1 -d " " | cut -f2 -d ">" > GRCh38_pop1.txt
Rscript ../../scripts/align_bowtie.R -g GRCh38_latest_rna -m 2 -e Cas13a -v GRCh38_pop1.txt

# align against cow transcriptome
Rscript ../../scripts/align_bowtie.R -g ARS-UCD1_rna -m 2 -e Cas13a

# calculate crRNA folding structures (RNAfold)
Rscript ../../scripts/score_RNAfold_crRNAs.R -e Cas13a

# compile scores and add extra annotations
Rscript compile_scores.R
