NC_045512v2.fa : wget from UCSC Genome Browser
http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/chromFa.tar.gz

wuhCor1 bowtie indices : bowtie-build NC_045512v2.fa wuhCor1

cov2_genes.txt : select fields from UCSC Genome Browser wuhCor1.ncbiGeneBGP

wuhCor1_RNAfold.txt : fold with RNAfold
RNAfold -i NC_045512v2.fa --outfile=wuhCor1_RNAfold.txt

wuhCor1.multiz119way.txt : download from UCSC Table Browser
https://genome-test.gi.ucsc.edu/cgi-bin/hgTables?hgsid=395041571_v0qlbJyAahP1qWEw9Yj5XfIz0lmK

wuhCor1.phyloP119way.txt : download from UCSC Table Browser
https://genome-test.gi.ucsc.edu/cgi-bin/hgTables?hgsid=395041571_v0qlbJyAahP1qWEw9Yj5XfIz0lmK&clade=virus&org=SARS-CoV-2&db=wuhCor1&hgta_group=compGeno&hgta_track=cons119way&hgta_table=phyloP119way&hgta_regionType=genome&position=NC_045512v2%3A1-29%2C903&hgta_outputType=maf&hgta_outFileName=

hg38.fa : wget from UCSC Genome Browser
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

hg38 bowtie indices : bowtie-build hg38.fa hg38

hg38_ncbiRefSeq.fa : download from UCSC Table Browser
clade: Mammal
genome: Human
assembly: Dec. 2013 (GRCh38/hg38)
group: Genes and Gene Predictions
track: NCBI RefSeq
table: RefSeq All (ncbiRefSeq)
output format: sequence
Sequence Retrieval Region Options:
- 5' UTR Exons
- CDS Exons
- 3' UTR Exons
- One FASTA record per gene
Sequence Formatting Options: All upper case

bosTau9.fa : wget from UCSC Genome Browser
http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/bosTau9.fa.gz

bosTau9 bowtie indices : bowtie-build bosTau9.fa bosTau9

bosTau9_ncbiRefSeq.fa : download from UCSC Table Browser
clade: Mammal
genome: Cow
assembly: Apr. 2018 (ARS-UCD1.2/bosTau9)
group: Genes and Gene Predictions
track: NCBI RefSeq
table: RefSeq All (ncbiRefSeq)
output format: sequence
Sequence Retrieval Region Options:
- 5' UTR Exons
- CDS Exons
- 3' UTR Exons
- One FASTA record per gene
Sequence Formatting Options: All upper case

EUA-Covid19-Template_0.docx : download from FDA
https://www.fda.gov/media/135658/download
- recommended list of organisms to be analyzed in silico

gisaid_cov2020_sequences.fasta : download from GISAID (2020/04/14 ; 5628 samples)
https://www.epicov.org/epi3/frontend#193bb6
- complete (>29,000bp)
- high coverage only

gisaid_cov2020_alignment.txt : generate pairwise alignment to wuhCor1 [Biostrings::pairwiseAlignment()]
Rscript ../scripts/align_gisaid_cov2020.R
sed -i 's/Hong Kong/HongKong/g' gisaid_cov2020_alignment.txt
sed -i 's/South Africa/SouthAfrica/g' gisaid_cov2020_alignment.txt
sed -i 's/South Korea/SouthKorea/g' gisaid_cov2020_alignment.txt
sed -i 's/Not found/NotFound/g' gisaid_cov2020_alignment.txt
sed -i 's/New Zealand/NewZealand/g' gisaid_cov2020_alignment.txt
sed -i 's/Lyon_06464 /Lyon_06464/' gisaid_cov2020_alignment.txt
sed -i 's/Czech Republic/CzechRepublic/g' gisaid_cov2020_alignment.txt
sed -i 's/Saudi Arabia/SaudiArabia/g' gisaid_cov2020_alignment.txt

# human coronaviruses (not SARS-CoV-2)

human_CoV/NC_004718v3.fa : SARS coronavirus, complete genome ; download from GenBank
https://www.ncbi.nlm.nih.gov/nuccore/NC_004718.3?report=fasta

human_CoV/NC_019843v3.fa : MERS coronavirus isolate, complete genome ; download from GenBank
https://www.ncbi.nlm.nih.gov/nuccore/NC_019843.3?report=fasta

human_CoV/NC_006213v1.fa : Human CoV OC43 strain ATCC VR-749 ; download from GenBank
https://www.ncbi.nlm.nih.gov/nuccore/NC_006213.1?report=fasta

human_CoV/NC_006577v2.fa : Human CoV HKU1 ; download from GenBank
https://www.ncbi.nlm.nih.gov/nuccore/NC_006577.2?report=fasta

human_CoV/NC_002645v1.fa : Human CoV 229E ; download from GenBank
https://www.ncbi.nlm.nih.gov/nuccore/NC_002645.1?report=fasta

human_CoV/NC_005831v2.fa : Human Coronavirus NL63 ; download from GenBank
https://www.ncbi.nlm.nih.gov/nuccore/NC_005831.2?report=fasta

cat human_CoV/*fa > human_CoV.fa

human_CoV bowtie indices: bowtie-build human_CoV.fa human_CoV

# human and cow transcriptomes

GRCh38_latest_rna.fna.gz : download from NCBI
ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz

GRCh38_latest_rna bowtie indices : bowtie-build GRCh38_latest_rna.fna GRCh38_latest_rna

GCF_002263795.1_ARS-UCD1.2_rna.fna.gz : download from NCBI
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_rna.fna.gz

ARS-UCD1_rna bowtie indices : bowtie-build GCF_002263795.1_ARS-UCD1.2_rna.fna ARS-UCD1_rna

# co-occurring pathogens

cat cross_reactive/*fa > cross_reactive.fa

cross_reactive bowtie indices : bowtie-build cross_reactive.fa cross_reactive

# clades from nextstrain

wuhCor1.nextstrainClade.txt : download from UCSC genome browser
clade: Viruses
genome: SARS-CoV-2
assembly: Jan. 2020/NC_045512.2
group: Variation and Repeats
track: Nextstrain Clades
table: nextstrainClade
region: genome
output format: selected fields from primary and related tables
selected fields: name, chromStarts, variants, sampleCount, samples

# viral genome structure

manfredonia_2020_lowShannon_highSHAPE.xlsx: download from Manfredonia (2020)
https://academic.oup.com/nar/article/48/22/12436/5961787#supplementary-data
Supplementary Table 2

MN985325v1.fa: SARS-CoV-2 isolate (USA-WA1)
lan_2020_structured.csv: personal correspondance
lan_2020_unstructured.csv: personal correspondance
bedtools getfasta -fi MN985325v1.fa -bed dms_map_unstructured.bed -fo dms_map_unstructured.fa
bowtie --norc -v 2 -S --un dms_map_unstructured_unmapped.fa -f wuhCor1 dms_map_unstructured.fa > dms_map_unstructured_mapped.sam 2> dms_map_unstructured_mapped.bowtiestats

huston_2021_structure.ct: download from Huston (2021)
https://github.com/pylelab/SARS-CoV-2_SHAPE_MaP_structure

icSHAPE_sun_2021.xlsx: download from Sun (2021)

influenza sequences downloaded from fludb.org
influenzaA: all subtypes, human host
influenzaB: human host

gisaid_omicron_variants_20211201.fasta
n = 309 viruses labeled "VOC Omicron GR/484A" from GISAID downloaded 1-dec-2021
