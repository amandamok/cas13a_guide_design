#!/bin/bash

# generate consensus sequences
for file in $(ls *fasta | cut -f1 -d".")
do
	echo "processing" ${file}
	usearch -cluster_fast ${file}.fasta -id 0.99 -centroids ${file}_centroids.fa
	muscle -in ${file}_centroids.fa -maxhours 3 -quiet -clwout ${file}.aln
	sed -i 's/MUSCLE/CLUSTAL/' ${file}.aln
	Rscript ../../scripts/generate_consensus.R -i ${file}.aln -n ${file} -o ${file}_consensus.fa
done

# aggregate segments
for num in $(seq 1 1 8)
do
	cat influenzaA_human_segment${num}_consensus.fa >> influenzaA_human_allSegments_consensus.fa
	cat influenzaB_human_segment${num}_consensus.fa >> influenzaB_human_allSegments_consensus.fa
	cat influenzaA_human_segment${num}_gene_consensus.fa >> influenzaA_human_allGenes_consensus.fa
	cat influenzaB_human_segment${num}_gene_consensus.fa >> influenzaB_human_allGenes_consensus.fa
done