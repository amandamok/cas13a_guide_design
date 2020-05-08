#!/bin/bash

### adenovirus: dsDNA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.gff.gz \
        -O adenovirus.gff.gz
gunzip adenovirus.gff.gz
grep -v ^"#" adenovirus.gff | cut -f1,4,5,7,9 | grep "ID=rna" | grep -v "tRNA-" | sed 's/;.*//g' | sed 's/ID=rna-//g' | \
        awk '{print $1, $2, $3, $5, $4}' OFS="\t" > adenovirus_rna.bed
for ID in $(cut -f1 adenovirus_rna.bed | sort | uniq)
do
        esearch -db nucleotide -query "${ID}" | efetch -format fasta >> adenovirus.fa
done
bedtools getfasta -s -fi adenovirus.fa -bed adenovirus_rna.bed -fo adenovirus_rna.fa

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

### respiratory syncytial virus A: ssRNA (-)
esearch -db nucleotide -query NC_038235.1 | efetch -format fasta > resp_syncytial_A.fa

### respiratory syncytial virus B: ssRNA (-)
esearch -db nucleotide -query NC_001781.1 | efetch -format fasta > resp_syncytial_B.fa

### rhinovirus: ssRNA (+)
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
