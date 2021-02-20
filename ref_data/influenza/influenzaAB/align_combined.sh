#!/bin/bash

for num in $(seq 1 1 8)
do
  file=$(echo "influenzaAB_human_segment"${num})
  echo "processing segment" ${file}
  cat ../influenzaA_human_segment${num}.fasta ../influenzaB_human_segment${num}.fasta \
      > ${file}.fa
  usearch -cluster_fast ${file}.fa -id 0.99 -centroids ${file}_centroids.fa
  muscle -in ${file}_centroids.fa -maxhours 3 -quiet -clwout ${file}.aln
  sed -i 's/MUSCLE/CLUSTAL/' ${file}.aln
  Rscript ../../../scripts/generate_consensus.R -i ${file}.aln \
      -n ${file} -o ${file}_consensus.fa
done
