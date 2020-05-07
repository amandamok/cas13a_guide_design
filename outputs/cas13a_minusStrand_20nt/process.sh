#!/bin/bash

rm *txt
rm *bowtiestats
rm *sam
rm *fa
rm *html

# break genome into windows of 20nt
Rscript ~/covid-19/scripts/break_genome.R -w 20 -e Cas13a -s -

# score sensitivity to SARS-CoV-2 genomes
Rscript ~/covid-19/scripts/score_SARS-CoV-2_sensitivity.R -w 20 -e Cas13a

# score specificity against other human coronaviruses
Rscript ~/covid-19/scripts/score_human_CoV_specificity.R -m 2

# calculate crRNA folding structures (RNAfold)
Rscript ~/covid-19/scripts/score_RNAfold_crRNAs.R -e Cas13a

# calculate viral target base-pairing propensity
Rscript ~/covid-19/scripts/score_RNAfold_target.R -w 20 -f 100

# off-target: align against hg38 (human)
Rscript ~/covid-19/scripts/align_bowtie.R -g GRCh38_latest_rna -m 2 -e Cas13a

# off-target: align against bosTau9 (cow)
Rscript ~/covid-19/scripts/align_bowtie.R -g ARS-UCD1_rna -m 2 -e Cas13a

# off-target: align against co-occuring pathogens
# Rscript ~/covid-19/scripts/align_bowtie.R -g cross_reactive

# compile results
Rscript compile_scores.R

# summarize results
