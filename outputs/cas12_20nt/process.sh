#!/bin/bash

rm *txt
rm *bowtiestats
rm *sam
rm *fa
rm *html

# break genome into windows of 20nt
Rscript ~/covid-19/scripts/break_genome.R -w 20 -e Cas12

# score sensitivity to SARS-CoV-2 genomes
Rscript ~/covid-19/scripts/score_SARS-CoV-2_sensitivity.R -w 20

# score specificity against other human coronaviruses
Rscript ~/covid-19/scripts/score_human_CoV_specificity.R -m 2 -e Cas12

# calculate crRNA folding structures (RNAfold)
Rscript ~/covid-19/scripts/score_RNAfold_crRNAs.R -e Cas12

# calculate viral target base-pairing propensity
Rscript ~/covid-19/scripts/score_RNAfold_target.R -w 20 -f 100

# off-target: align against hg38 (human)
Rscript ~/covid-19/scripts/align_bowtie.R -g hg38 -m 2 -e Cas12

# off-target: align against bosTau9 (cow)
Rscript ~/covid-19/scripts/align_bowtie.R -g bosTau9 -m 2 -e Cas12

# off-target: align against co-occuring pathogens
# Rscript ~/covid-19/scripts/align_bowtie.R -g cross_reactive

# compile results
Rscript compile_scores.R

# summarize results
Rscript -e 'rmarkdown::render("cas12_results_summary_test20.Rmd", quiet=T)' 
