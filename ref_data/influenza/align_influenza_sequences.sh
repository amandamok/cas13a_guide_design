#!/bin/bash

# generate consensus sequences
for file in $(ls *fasta | cut -f1 -d".")
do
	echo "processing" ${file}
	usearch -cluster_fast ${file}.fasta -id 0.99 -centroids ${file}_centroids.fa
	muscle -in ${file}_centroids.fa -maxhours 3 -quiet -clwout ${file}.aln
	sed -i 's/MUSCLE/CLUSTAL/' ${file}.aln
	Rscript ../../scripts/generate_consensus.R -i ${file}.aln -t clustal -n ${file} -o ${file}_consensus.fa
done

# aggregate segments
for num in $(seq 1 1 8)
do
	cat influenzaA_human_segment${num}_consensus.fa >> influenzaA_human_allSegments_consensus.fa
	cat influenzaB_human_segment${num}_consensus.fa >> influenzaB_human_allSegments_consensus.fa
	cat influenzaA_human_segment${num}_gene_consensus.fa >> influenzaA_human_allGenes_consensus.fa
	cat influenzaB_human_segment${num}_gene_consensus.fa >> influenzaB_human_allGenes_consensus.fa
done

# influenzaA_H1N1_human_2015-2020
cd influenzaA_H1N1_human_2015-2020
for file in $(ls *fasta | cut -f1 -d".")
do
    Rscript ../../../scripts/generate_consensus.R -i ${file}.fasta -t fasta \
        -n ${file} -o ${file}_consensus.fa
    cat ${file}_consensus.fa >> influenzaA_H1N1_human_2015-2020_allSegments_consensus.fa
done

# influenzaA_H3N2_human_2015-2020
cd ../influenzaA_H3N2_human_2015-2020
for file in $(ls *fasta | cut -f1 -d".")
do
    Rscript ../../../scripts/generate_consensus.R -i ${file}.fasta -t fasta \
        -n ${file} -o ${file}_consensus.fa
    cat ${file}_consensus.fa >> influenzaA_H3N2_human_2015-2020_allSegments_consensus.fa
done

# influenzaB_human_2015-2021
cd ../influenzaB_human_2015-2021
for file in $(ls *fasta | cut -f1 -d".")
do
    Rscript ../../../scripts/generate_consensus.R -i ${file}.fasta -t fasta \
        -n ${file} -o ${file}_consensus.fa
    cat ${file}_consensus.fa >> influenzaB_human_2015-2021_allSegments_consensus.fa
done
