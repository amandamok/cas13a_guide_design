#!/bin/bash

rm *txt
rm *bowtiestats
rm *sam
rm *fa
rm *html
rm *pdf

# create amplicon fasta file
echo CCAAGTAATTGAAAAGACACTCCTCCACTTATCC > amplicon.fa

# break genome into windows of 20nt
Rscript ~/cas13a_guide_design/scripts/break_genome.R -i amplicon.fa -w 20 -e Cas13a

# score crRNA folding structures (RNAfold)
Rscript ~/cas13a_guide_design/scripts/score_RNAfold_crRNAs.R -e Cas13a

# off-target: align against hg38 (human)
Rscript ~/cas13a_guide_design/scripts/align_bowtie.R -g GRCh38_latest_rna -m 2 -e Cas13a

# off-target: align against SARS-CoV-2
Rscript ~/cas13a_guide_design/scripts/align_bowtie.R -g wuhCor1 -m 2 -e Cas13a

# compile results
Rscript compile_scores.R
