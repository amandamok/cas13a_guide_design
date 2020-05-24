#!/bin/bash

### hMPV: ssRNA (-)
esearch -db nucleotide -query NC_039199.1 | efetch -format fasta > hMPV.fa

### parainfluenza 1: ssRNA (-)
esearch -db nucleotide -query NC_003461.1 | efetch -format fasta > parainfluenza_1.fa

### parainfluenza 2: ssRNA (-)
esearch -db nucleotide -query NC_003443.1 | efetch -format fasta > parainfluenza_2.fa

### parainfluenza 3: ssRNA (-)
esearch -db nucleotide -query NC_001796.2 | efetch -format fasta > parainfluenza_3.fa

### parainfluenza 4: ssRNA (-)
esearch -db nucleotide -query NC_021928.1 | efetch -format fasta > parainfluenza_4.fa

### influenza A: ssRNA (-)
esearch -db nucleotide -query NC_007373.1 | efetch -format fasta > influenza_A.fa
esearch -db nucleotide -query NC_007372.1 | efetch -format fasta >> influenza_A.fa
esearch -db nucleotide -query NC_007371.1 | efetch -format fasta >> influenza_A.fa
esearch -db nucleotide -query NC_007366.1 | efetch -format fasta >> influenza_A.fa
esearch -db nucleotide -query NC_007369.1 | efetch -format fasta >> influenza_A.fa
esearch -db nucleotide -query NC_007368.1 | efetch -format fasta >> influenza_A.fa
esearch -db nucleotide -query NC_007367.1 | efetch -format fasta >> influenza_A.fa
esearch -db nucleotide -query NC_007370.1 | efetch -format fasta >> influenza_A.fa

### influenza B: ssRNA (-)
esearch -db nucleotide -query NC_002204.1 | efetch -format fasta > influenza_B.fa
esearch -db nucleotide -query NC_002205.1 | efetch -format fasta >> influenza_B.fa
esearch -db nucleotide -query NC_002206.1 | efetch -format fasta >> influenza_B.fa
esearch -db nucleotide -query NC_002207.1 | efetch -format fasta >> influenza_B.fa
esearch -db nucleotide -query NC_002208.1 | efetch -format fasta >> influenza_B.fa
esearch -db nucleotide -query NC_002209.1 | efetch -format fasta >> influenza_B.fa
esearch -db nucleotide -query NC_002210.1 | efetch -format fasta >> influenza_B.fa
esearch -db nucleotide -query NC_002211.1 | efetch -format fasta >> influenza_B.fa

### influenza C: ssRNA (-)
esearch -db nucleotide -query NC_006307.2 | efetch -format fasta > influenza_C.fa
esearch -db nucleotide -query NC_006308.2 | efetch -format fasta >> influenza_C.fa
esearch -db nucleotide -query NC_006309.2 | efetch -format fasta >> influenza_C.fa
esearch -db nucleotide -query NC_006310.2 | efetch -format fasta >> influenza_C.fa
esearch -db nucleotide -query NC_006311.1 | efetch -format fasta >> influenza_C.fa
esearch -db nucleotide -query NC_006312.2 | efetch -format fasta >> influenza_C.fa
esearch -db nucleotide -query NC_006306.2 | efetch -format fasta >> influenza_C.fa

### enterovirus: ssRNA (+)
esearch -db nucleotide -query NC_038308.1 | efetch -format fasta > enterovirus.fa

### respiratory syncytial virus A (human orthopneumovirus): ssRNA (-)
esearch -db nucleotide -query NC_038235.1 | efetch -format fasta > resp_syncytial_A.fa

### respiratory syncytial virus B (human orthopneumovirus): ssRNA (-)
esearch -db nucleotide -query NC_001781.1 | efetch -format fasta > resp_syncytial_B.fa

### rhinovirus C: ssRNA (+)
esearch -db nucleotide -query NC_009996.1 | efetch -format fasta > rhinovirus.fa

### parechovirus: ssRNA (+)
esearch -db nucleotide -query NC_001897.1 | efetch -format fasta > parechovirus.fa

### candida albicans: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.gff.gz \
	-O c_albicans.gff.gz
gunzip c_albicans.gff.gz
grep -v ^"#" c_albicans.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
	awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_albicans_rna.bed
for ID in $(cut -f1 c_albicans_rna.bed | sort | uniq)
do
	esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_albicans.fa
done
bedtools getfasta -s -fi c_albicans.fa -bed c_albicans_rna.bed -fo c_albicans_rna.fa

### corynebacterium diphtheriae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/457/455/GCF_001457455.1_NCTC11397/GCF_001457455.1_NCTC11397_genomic.gff.gz \
        -O c_diphtheriae.gff.gz
gunzip c_diphtheriae.gff.gz
grep -v ^"#" c_diphtheriae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_diphtheriae_rna.bed
for ID in $(cut -f1 c_diphtheriae_rna.bed | sort | uniq)
do      
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_diphtheriae.fa
done
bedtools getfasta -s -fi c_diphtheriae.fa -bed c_diphtheriae_rna.bed -fo c_diphtheriae_rna.fa

### legionella anisa: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/467/525/GCF_001467525.1_ASM146752v1/GCF_001467525.1_ASM146752v1_genomic.gff.gz \
        -O l_anisa.gff.gz
gunzip l_anisa.gff.gz
grep -v ^"#" l_anisa.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > l_anisa_rna.bed
for ID in $(cut -f1 l_anisa_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> l_anisa.fa
done
bedtools getfasta -s -fi l_anisa.fa -bed l_anisa_rna.bed -fo l_anisa_rna.fa

### legionella pneumophila: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/485/GCF_000008485.1_ASM848v1/GCF_000008485.1_ASM848v1_genomic.gff.gz \
        -O l_pneumophila.gff.gz
gunzip l_pneumophila.gff.gz
grep -v ^"#" l_pneumophila.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > l_pneumophila_rna.bed
for ID in $(cut -f1 l_pneumophila_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> l_pneumophila.fa
done
bedtools getfasta -s -fi l_pneumophila.fa -bed l_pneumophila_rna.bed -fo l_pneumophila_rna.fa

### bacillus anthracis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/845/GCF_000007845.1_ASM784v1/GCF_000007845.1_ASM784v1_genomic.gff.gz \
        -O b_anthracis.gff.gz
gunzip b_anthracis.gff.gz
grep -v ^"#" b_anthracis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > b_anthracis_rna.bed
for ID in $(cut -f1 b_anthracis_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> b_anthracis.fa
done
bedtools getfasta -s -fi b_anthracis.fa -bed b_anthracis_rna.bed -fo b_anthracis_rna.fa

### moraxella catarrhalis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/092/265/GCF_000092265.1_ASM9226v1/GCF_000092265.1_ASM9226v1_genomic.gff.gz \
        -O m_catarrhalis.gff.gz
gunzip m_catarrhalis.gff.gz
grep -v ^"#" m_catarrhalis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > m_catarrhalis_rna.bed
for ID in $(cut -f1 m_catarrhalis_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> m_catarrhalis.fa
done
bedtools getfasta -s -fi m_catarrhalis.fa -bed m_catarrhalis_rna.bed -fo m_catarrhalis_rna.fa

### neisseria elongata: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/818/035/GCF_000818035.1_ASM81803v1/GCF_000818035.1_ASM81803v1_genomic.gff.gz \
        -O n_elongata.gff.gz
gunzip n_elongata.gff.gz
grep -v ^"#" n_elongata.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > n_elongata_rna.bed
for ID in $(cut -f1 n_elongata_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> n_elongata.fa
done
bedtools getfasta -s -fi n_elongata.fa -bed n_elongata_rna.bed -fo n_elongata_rna.fa

### neisseria meningitidis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/805/GCF_000008805.1_ASM880v1/GCF_000008805.1_ASM880v1_genomic.gff.gz \
        -O n_meningitidis.gff.gz
gunzip n_meningitidis.gff.gz
grep -v ^"#" n_meningitidis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > n_meningitidis_rna.bed
for ID in $(cut -f1 n_meningitidis_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> n_meningitidis.fa
done
bedtools getfasta -s -fi n_meningitidis.fa -bed n_meningitidis_rna.bed -fo n_meningitidis_rna.fa

### pseudomonas aeruginosa: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz \
        -O p_aeruginosa.gff.gz
gunzip p_aeruginosa.gff.gz
grep -v ^"#" p_aeruginosa.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > p_aeruginosa_rna.bed
for ID in $(cut -f1 p_aeruginosa_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> p_aeruginosa.fa
done
bedtools getfasta -s -fi p_aeruginosa.fa -bed p_aeruginosa_rna.bed -fo p_aeruginosa_rna.fa

### staphylococcus epidermidis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/645/GCF_000007645.1_ASM764v1/GCF_000007645.1_ASM764v1_genomic.gff.gz \
        -O s_epidermidis.gff.gz
gunzip s_epidermidis.gff.gz
grep -v ^"#" s_epidermidis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_epidermidis_rna.bed
for ID in $(cut -f1 s_epidermidis_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_epidermidis.fa
done
bedtools getfasta -s -fi s_epidermidis.fa -bed s_epidermidis_rna.bed -fo s_epidermidis_rna.fa

### streptococcus salivarius: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/785/515/GCF_000785515.1_ASM78551v1/GCF_000785515.1_ASM78551v1_genomic.gff.gz \
        -O s_salivarius.gff.gz
gunzip s_salivarius.gff.gz
grep -v ^"#" s_salivarius.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_salivarius_rna.bed
for ID in $(cut -f1 s_salivarius_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_salivarius.fa
done
bedtools getfasta -s -fi s_salivarius.fa -bed s_salivarius_rna.bed -fo s_salivarius_rna.fa

### leptospira santarosai: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/313/175/GCF_000313175.2_ASM31317v2/GCF_000313175.2_ASM31317v2_genomic.gff.gz \
        -O l_santarosai.gff.gz
gunzip l_santarosai.gff.gz
grep -v ^"#" l_santarosai.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > l_santarosai_rna.bed
for ID in $(cut -f1 l_santarosai_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> l_santarosai.fa
done
bedtools getfasta -s -fi l_santarosai.fa -bed l_santarosai_rna.bed -fo l_santarosai_rna.fa

### chlamydia pneumoniae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/745/GCF_000008745.1_ASM874v1/GCF_000008745.1_ASM874v1_genomic.gff.gz \
        -O c_pneumoniae.gff.gz
gunzip c_pneumoniae.gff.gz
grep -v ^"#" c_pneumoniae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_pneumoniae_rna.bed
for ID in $(cut -f1 c_pneumoniae_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_pneumoniae.fa
done
bedtools getfasta -s -fi c_pneumoniae.fa -bed c_pneumoniae_rna.bed -fo c_pneumoniae_rna.fa

### chlamydia psittaci: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/204/255/GCF_000204255.1_ASM20425v1/GCF_000204255.1_ASM20425v1_genomic.gff.gz \
        -O c_psittaci.gff.gz
gunzip c_psittaci.gff.gz
grep -v ^"#" c_psittaci.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_psittaci_rna.bed
for ID in $(cut -f1 c_psittaci_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_psittaci.fa
done
bedtools getfasta -s -fi c_psittaci.fa -bed c_psittaci_rna.bed -fo c_psittaci_rna.fa

### coxiella burnetii: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/765/GCF_000007765.2_ASM776v2/GCF_000007765.2_ASM776v2_genomic.gff.gz \
        -O c_burnetii.gff.gz
gunzip c_burnetii.gff.gz
grep -v ^"#" c_burnetii.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_burnetii_rna.bed
for ID in $(cut -f1 c_burnetii_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_burnetii.fa
done
bedtools getfasta -s -fi c_burnetii.fa -bed c_burnetii_rna.bed -fo c_burnetii_rna.fa

### stapylococcus aureus: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz \
        -O s_aureus.gff.gz
gunzip s_aureus.gff.gz
grep -v ^"#" s_aureus.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_aureus_rna.bed
for ID in $(cut -f1 s_aureus_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_aureus.fa
done
bedtools getfasta -s -fi s_aureus.fa -bed s_aureus_rna.bed -fo s_aureus_rna.fa

### haemophilus influenzae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/305/GCF_000027305.1_ASM2730v1/GCF_000027305.1_ASM2730v1_genomic.gff.gz \
        -O h_influenzae.gff.gz
gunzip h_influenzae.gff.gz
grep -v ^"#" h_influenzae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > h_influenzae_rna.bed
for ID in $(cut -f1 h_influenzae_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> h_influenzae.fa
done
bedtools getfasta -s -fi h_influenzae.fa -bed h_influenzae_rna.bed -fo h_influenzae_rna.fa

### mycobacterium tuberculosis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz \
        -O m_tuberculosis.gff.gz
gunzip m_tuberculosis.gff.gz
grep -v ^"#" m_tuberculosis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > m_tuberculosis_rna.bed
for ID in $(cut -f1 m_tuberculosis_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> m_tuberculosis.fa
done
bedtools getfasta -s -fi m_tuberculosis.fa -bed m_tuberculosis_rna.bed -fo m_tuberculosis_rna.fa

### streptococcus pneumoniae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/045/GCF_000007045.1_ASM704v1/GCF_000007045.1_ASM704v1_genomic.gff.gz \
        -O s_pneumoniae.gff.gz
gunzip s_pneumoniae.gff.gz
grep -v ^"#" s_pneumoniae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_pneumoniae_rna.bed
for ID in $(cut -f1 s_pneumoniae_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_pneumoniae.fa
done
bedtools getfasta -s -fi s_pneumoniae.fa -bed s_pneumoniae_rna.bed -fo s_pneumoniae_rna.fa

### streptococcus pyogenes: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/785/GCF_000006785.2_ASM678v2/GCF_000006785.2_ASM678v2_genomic.gff.gz \
        -O s_pyogenes.gff.gz
gunzip s_pyogenes.gff.gz
grep -v ^"#" s_pyogenes.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_pyogenes_rna.bed
for ID in $(cut -f1 s_pyogenes_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_pyogenes.fa
done
bedtools getfasta -s -fi s_pyogenes.fa -bed s_pyogenes_rna.bed -fo s_pyogenes_rna.fa

### bordetella pertussis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/715/GCF_000195715.1_ASM19571v1/GCF_000195715.1_ASM19571v1_genomic.gff.gz \
        -O b_pertussis.gff.gz
gunzip b_pertussis.gff.gz
grep -v ^"#" b_pertussis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > b_pertussis_rna.bed
for ID in $(cut -f1 b_pertussis_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> b_pertussis.fa
done
bedtools getfasta -s -fi b_pertussis.fa -bed b_pertussis_rna.bed -fo b_pertussis_rna.fa

### mycoplasma pneumoniae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/345/GCF_000027345.1_ASM2734v1/GCF_000027345.1_ASM2734v1_genomic.gff.gz \
        -O m_pneumoniae.gff.gz
gunzip m_pneumoniae.gff.gz
grep -v ^"#" m_pneumoniae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > m_pneumoniae_rna.bed
for ID in $(cut -f1 m_pneumoniae_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> m_pneumoniae.fa
done
bedtools getfasta -s -fi m_pneumoniae.fa -bed m_pneumoniae_rna.bed -fo m_pneumoniae_rna.fa

### pneumocystis jirovecii: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/477/535/GCF_001477535.1_Pneu_jiro_RU7_V2/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.gff.gz \
        -O p_jirovecii.gff.gz
gunzip p_jirovecii.gff.gz
grep -v ^"#" p_jirovecii.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > p_jirovecii_rna.bed
for ID in $(cut -f1 p_jirovecii_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> p_jirovecii.fa
done
bedtools getfasta -s -fi p_jirovecii.fa -bed p_jirovecii_rna.bed -fo p_jirovecii_rna.fa

### streptococcus intermedius: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/463/355/GCF_000463355.1_ASM46335v1/GCF_000463355.1_ASM46335v1_genomic.gff.gz \
        -O s_intermedius.gff.gz
gunzip s_intermedius.gff.gz
grep -v ^"#" s_intermedius.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_intermedius_rna.bed
for ID in $(cut -f1 s_intermedius_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_intermedius.fa
done
bedtools getfasta -s -fi s_intermedius.fa -bed s_intermedius_rna.bed -fo s_intermedius_rna.fa

### streptococcus constellatus: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/463/425/GCF_000463425.1_ASM46342v1/GCF_000463425.1_ASM46342v1_genomic.gff.gz \
      -O s_constellatus.gff.gz
gunzip s_constellatus.gff.gz
grep -v ^"#" s_constellatus.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_constellatus_rna.bed
for ID in $(cut -f1 s_constellatus_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_constellatus.fa
done
bedtools getfasta -s -fi s_constellatus.fa -bed s_constellatus_rna.bed -fo s_constellatus_rna.fa

### streptococcus anginosus: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/463/505/GCF_000463505.1_ASM46350v1/GCF_000463505.1_ASM46350v1_genomic.gff.gz \
      -O s_anginosus.gff.gz
gunzip s_anginosus.gff.gz
grep -v ^"#" s_anginosus.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_anginosus_rna.bed
for ID in $(cut -f1 s_anginosus_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_anginosus.fa
done
bedtools getfasta -s -fi s_anginosus.fa -bed s_anginosus_rna.bed -fo s_anginosus_rna.fa

### stenotrophomonas maltophilia: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/072/485/GCF_000072485.1_ASM7248v1/GCF_000072485.1_ASM7248v1_genomic.gff.gz \
      -O s_maltophilia.gff.gz
gunzip s_maltophilia.gff.gz
grep -v ^"#" s_maltophilia.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_maltophilia_rna.bed
for ID in $(cut -f1 s_maltophilia_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_maltophilia.fa
done
bedtools getfasta -s -fi s_maltophilia.fa -bed s_maltophilia_rna.bed -fo s_maltophilia_rna.fa

### serratia marcescens: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/513/215/GCF_000513215.1_DB11/GCF_000513215.1_DB11_genomic.gff.gz \
      -O s_marcescens.gff.gz
gunzip s_marcescens.gff.gz
grep -v ^"#" s_marcescens.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > s_marcescens_rna.bed
for ID in $(cut -f1 s_marcescens_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> s_marcescens.fa
done
bedtools getfasta -s -fi s_marcescens.fa -bed s_marcescens_rna.bed -fo s_marcescens_rna.fa

### proteus mirabilis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/069/965/GCF_000069965.1_ASM6996v1/GCF_000069965.1_ASM6996v1_genomic.gff.gz \
      -O p_mirabilis.gff.gz
gunzip p_mirabilis.gff.gz
grep -v ^"#" p_mirabilis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > p_mirabilis_rna.bed
for ID in $(cut -f1 p_mirabilis_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> p_mirabilis.fa
done
bedtools getfasta -s -fi p_mirabilis.fa -bed p_mirabilis_rna.bed -fo p_mirabilis_rna.fa

### pasteurella multocida: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/754/275/GCF_000754275.1_ASM75427v1/GCF_000754275.1_ASM75427v1_genomic.gff.gz \
      -O p_multocida.gff.gz
gunzip p_multocida.gff.gz
grep -v ^"#" p_multocida.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > p_multocida_rna.bed
for ID in $(cut -f1 p_multocida_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> p_multocida.fa
done
bedtools getfasta -s -fi p_multocida.fa -bed p_multocida_rna.bed -fo p_multocida_rna.fa

### morganella morganii: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/286/435/GCF_000286435.2_ASM28643v2/GCF_000286435.2_ASM28643v2_genomic.gff.gz \
      -O m_morganii.gff.gz
gunzip m_morganii.gff.gz
grep -v ^"#" m_morganii.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > m_morganii_rna.bed
for ID in $(cut -f1 m_morganii_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> m_morganii.fa
done
bedtools getfasta -s -fi m_morganii.fa -bed m_morganii_rna.bed -fo m_morganii_rna.fa

### klebsiella pneumoniae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.gff.gz \
      -O k_pneumoniae.gff.gz
gunzip k_pneumoniae.gff.gz
grep -v ^"#" k_pneumoniae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > k_pneumoniae_rna.bed
for ID in $(cut -f1 k_pneumoniae_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> k_pneumoniae.fa
done
bedtools getfasta -s -fi k_pneumoniae.fa -bed k_pneumoniae_rna.bed -fo k_pneumoniae_rna.fa

### klebsiella oxytoca: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/022/195/GCF_001022195.1_ASM102219v1/GCF_001022195.1_ASM102219v1_genomic.gff.gz \
      -O k_oxytoca.gff.gz
gunzip k_oxytoca.gff.gz
grep -v ^"#" k_oxytoca.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > k_oxytoca_rna.bed
for ID in $(cut -f1 k_oxytoca_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> k_oxytoca.fa
done
bedtools getfasta -s -fi k_oxytoca.fa -bed k_oxytoca_rna.bed -fo k_oxytoca_rna.fa

### klebsiella aerogenes: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/215/745/GCF_000215745.1_ASM21574v1/GCF_000215745.1_ASM21574v1_genomic.gff.gz \
      -O k_aerogenes.gff.gz
gunzip k_aerogenes.gff.gz
grep -v ^"#" k_aerogenes.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > k_aerogenes_rna.bed
for ID in $(cut -f1 k_aerogenes_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> k_aerogenes.fa
done
bedtools getfasta -s -fi k_aerogenes.fa -bed k_aerogenes_rna.bed -fo k_aerogenes_rna.fa

### fusobacterium nucleatum: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/325/GCF_000007325.1_ASM732v1/GCF_000007325.1_ASM732v1_genomic.gff.gz \
      -O f_nucleatum.gff.gz
gunzip f_nucleatum.gff.gz
grep -v ^"#" f_nucleatum.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > f_nucleatum_rna.bed
for ID in $(cut -f1 f_nucleatum_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> f_nucleatum.fa
done
bedtools getfasta -s -fi f_nucleatum.fa -bed f_nucleatum_rna.bed -fo f_nucleatum_rna.fa

### fusobacterium necrophorum: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/262/225/GCF_000262225.1_ASM26222v1/GCF_000262225.1_ASM26222v1_genomic.gff.gz \
      -O f_necrophorum.gff.gz
gunzip f_necrophorum.gff.gz
grep -v ^"#" f_necrophorum.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > f_necrophorum_rna.bed
for ID in $(cut -f1 f_necrophorum_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> f_necrophorum.fa
done
bedtools getfasta -s -fi f_necrophorum.fa -bed f_necrophorum_rna.bed -fo f_necrophorum_rna.fa

### francisella tularensis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/985/GCF_000008985.1_ASM898v1/GCF_000008985.1_ASM898v1_genomic.gff.gz \
      -O f_tularensis.gff.gz
gunzip f_tularensis.gff.gz
grep -v ^"#" f_tularensis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > f_tularensis_rna.bed
for ID in $(cut -f1 f_tularensis_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> f_tularensis.fa
done
bedtools getfasta -s -fi f_tularensis.fa -bed f_tularensis_rna.bed -fo f_tularensis_rna.fa

### escherichia coli: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz \
      -O e_coli.gff.gz
gunzip e_coli.gff.gz
grep -v ^"#" e_coli.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > e_coli_rna.bed
for ID in $(cut -f1 e_coli_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> e_coli.fa
done
bedtools getfasta -s -fi e_coli.fa -bed e_coli_rna.bed -fo e_coli_rna.fa

### enterobacter cloacae: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/565/GCF_000025565.1_ASM2556v1/GCF_000025565.1_ASM2556v1_genomic.gff.gz \
      -O e_cloacae.gff.gz
gunzip e_cloacae.gff.gz
grep -v ^"#" e_cloacae.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > e_cloacae_rna.bed
for ID in $(cut -f1 e_cloacae_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> e_cloacae.fa
done
bedtools getfasta -s -fi e_cloacae.fa -bed e_cloacae_rna.bed -fo e_cloacae_rna.fa

### citrobacter koseri: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/045/GCF_000018045.1_ASM1804v1/GCF_000018045.1_ASM1804v1_genomic.gff.gz \
      -O c_koseri.gff.gz
gunzip c_koseri.gff.gz
grep -v ^"#" c_koseri.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_koseri_rna.bed
for ID in $(cut -f1 c_koseri_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_koseri.fa
done
bedtools getfasta -s -fi c_koseri.fa -bed c_koseri_rna.bed -fo c_koseri_rna.fa

### citrobacter freundii: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/648/515/GCF_000648515.1_ASM64851v1/GCF_000648515.1_ASM64851v1_genomic.gff.gz \
      -O c_freundii.gff.gz
gunzip c_freundii.gff.gz
grep -v ^"#" c_freundii.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > c_freundii_rna.bed
for ID in $(cut -f1 c_freundii_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> c_freundii.fa
done
bedtools getfasta -s -fi c_freundii.fa -bed c_freundii_rna.bed -fo c_freundii_rna.fa

### burkholderia pseudomallei: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/545/GCF_000011545.1_ASM1154v1/GCF_000011545.1_ASM1154v1_genomic.gff.gz \
      -O b_pseudomallei.gff.gz
gunzip b_pseudomallei.gff.gz
grep -v ^"#" b_pseudomallei.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > b_pseudomallei_rna.bed
for ID in $(cut -f1 b_pseudomallei_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> b_pseudomallei.fa
done
bedtools getfasta -s -fi b_pseudomallei.fa -bed b_pseudomallei_rna.bed -fo b_pseudomallei_rna.fa

### bacteroides fragilis: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/925/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.gff.gz \
      -O b_fragilis.gff.gz
gunzip b_fragilis.gff.gz
grep -v ^"#" b_fragilis.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > b_fragilis_rna.bed
for ID in $(cut -f1 b_fragilis_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> b_fragilis.fa
done
bedtools getfasta -s -fi b_fragilis.fa -bed b_fragilis_rna.bed -fo b_fragilis_rna.fa

### acinetobacter baumannii: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/746/645/GCF_000746645.1_ASM74664v1/GCF_000746645.1_ASM74664v1_genomic.gff.gz \
      -O a_baumannii.gff.gz
gunzip a_baumannii.gff.gz
grep -v ^"#" a_baumannii.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > a_baumannii_rna.bed
for ID in $(cut -f1 a_baumannii_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> a_baumannii.fa
done
bedtools getfasta -s -fi a_baumannii.fa -bed a_baumannii_rna.bed -fo a_baumannii_rna.fa

### aspergillus niger: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/855/GCF_000002855.3_ASM285v2/GCF_000002855.3_ASM285v2_genomic.gff.gz \
      -O a_niger.gff.gz
gunzip a_niger.gff.gz
grep -v ^"#" a_niger.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > a_niger_rna.bed
for ID in $(cut -f1 a_niger_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> a_niger.fa
done
bedtools getfasta -s -fi a_niger.fa -bed a_niger_rna.bed -fo a_niger_rna.fa

### aspergillus flavus: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/275/GCF_000006275.2_JCVI-afl1-v2.0/GCF_000006275.2_JCVI-afl1-v2.0_genomic.gff.gz \
      -O a_flavus.gff.gz
gunzip a_flavus.gff.gz
grep -v ^"#" a_flavus.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
      awk '{print $1, $2, $3, $5, $4}' OFS="\t" > a_flavus_rna.bed
for ID in $(cut -f1 a_flavus_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> a_flavus.fa
done
bedtools getfasta -s -fi a_flavus.fa -bed a_flavus_rna.bed -fo a_flavus_rna.fa

### hCMV (human betaherpesvirus 5): dsDNA
# download gene table from https://www.ncbi.nlm.nih.gov/gene/?term=txid10359[Organism:noexp] > hCMV.gene
cut -f12,13,14,15,6 hCMV.gene | tail -n +2 | head -n -1 | awk '{print $2,$3,$4,$5,$1}' OFS="\t" \
	| sed 's/plus/+/g' | sed 's/minus/-/g' > hCMV_rna.bed
for ID in $(cut -f1 hCMV_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> hCMV.fa
done
bedtools getfasta -s -fi hCMV.fa -bed hCMV_rna.bed -fo hCMV_rna.fa

### human mastadenovirus D: dsDNA
# download gene table from https://www.ncbi.nlm.nih.gov/gene?LinkName=genome_gene&from_uid=10278 > h_mastadenovirus_D.gene
cut -f12,13,14,15,6 h_mastadenovirus_D.gene | tail -n +2 | head -n -1 | awk '{print $2,$3,$4,$5,$1}' OFS="\t" \
	| sed 's/plus/+/g' | sed 's/minus/-/g' > h_mastadenovirus_D_rna.bed
for ID in $(cut -f1 h_mastadenovirus_D_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> h_mastadenovirus_D.fa
done
bedtools getfasta -s -fi h_mastadenovirus_D.fa -bed h_mastadenovirus_D_rna.bed -fo h_mastadenovirus_D_rna.fa

### human mastadenovirus C: dsDNA
# download gene table from https://www.ncbi.nlm.nih.gov/gene?LinkName=genome_gene&from_uid=10271 > h_mastadenovirus_C.gene
cut -f12,13,14,15,6 h_mastadenovirus_C.gene | tail -n +2 | head -n -1 | awk '{print $2,$3,$4,$5,$1}' OFS="\t" \
	| sed 's/plus/+/g' | sed 's/minus/-/g' > h_mastadenovirus_C_rna.bed
for ID in $(cut -f1 h_mastadenovirus_C_rna.bed | sort | uniq)
do
      esearch -db nucleotide -query "${ID}" | efetch -format fasta >> h_mastadenovirus_C.fa
done
bedtools getfasta -s -fi h_mastadenovirus_C.fa -bed h_mastadenovirus_C_rna.bed -fo h_mastadenovirus_C_rna.fa
