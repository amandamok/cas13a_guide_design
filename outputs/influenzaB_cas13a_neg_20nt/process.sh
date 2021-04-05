#!/bin/bash

rm *txt
rm *bowtiestats
rm *sam
rm *fa

script_dir="../../scripts"
ref_dir="../../ref_data/influenza"

# break genome segments into windows of 20nt
Rscript ${script_dir}/break_genome.R -i ${ref_dir}/influenzaB_human_allSegments_consensus.fa \
    -w 20 -e Cas13a -g - -s -

# score sensitivity to reference genome segments
Rscript ${script_dir}/score_influenza_sensitivity.R -t influenzaB -m 1 -w 20

# score crRNA folding structures (RNAfold)
Rscript ${script_dir}/score_RNAfold_crRNAs.R -e Cas13a

# off-target: align against hg38 (human)
Rscript ${script_dir}/align_bowtie.R -g GRCh38_latest_rna -m 2 -e Cas13a

# off-target: align against bosTau9 (cow)
Rscript ${script_dir}/align_bowtie.R -g ARS-UCD1_rna -m 2 -e Cas13a

# compile results
Rscript compile_scores.R