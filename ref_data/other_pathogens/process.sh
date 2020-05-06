#!/bin/bash

# paste viral accession ids
# copy over .csv and .fa files over to directory

cut -f2 -d, *csv | sort | uniq | grep -v Accession >> accession_ids.txt
sed -i 's/"//g' accession_ids.txt

for ID in $(cat accession_ids.txt)
do
    echo $ID
    esearch -db nucleotide -query "$ID" | efetch -format fasta > $ID.fa
done


