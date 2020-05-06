#!/bin/bash

# SRR11140750: paired, clinical (nasal?) swab from Madison, WI
Rscript ~/covid-19/scripts/align_viral_abundance.R -a SRR11140750 -p

# SRR11454606: paired, throat swab
Rscript ~/covid-19/scripts/align_viral_abundance.R -a SRR11454606 -p

# SRR11454609: single, throat swab
Rscript ~/covid-19/scripts/align_viral_abundance.R -a SRR11454609

# SRR11454610: single, throat swab
Rscript ~/covid-19/scripts/align_viral_abundance.R -a SRR11454610

# SRR11454611: single, throat swab
Rscript ~/covid-19/scripts/align_viral_abundance.R -a SRR11454611

# concatenate SRR11454606, SRR11454609, SRR11454610, SRR11454611
cat SRR11454606_mapped.sam SRR11454609_mapped.sam SRR11454610_mapped.sam SRR11454611_mapped.sam > PRJNA616446_mapped.sam
Rscript ~/covid-19/scripts/align_viral_abundance.R -a PRJNA616446

# compute coverage with bedtools genomecov
samtools view -b SRR11454606_mapped.sam > SRR11454606_mapped.bam
samtools view -b SRR11454609_mapped.sam > SRR11454609_mapped.bam
samtools view -b SRR11454610_mapped.sam > SRR11454610_mapped.bam
samtools view -b SRR11454611_mapped.sam > SRR11454611_mapped.bam
samtools cat SRR11454606_mapped.bam SRR11454609_mapped.bam SRR11454610_mapped.bam SRR11454611_mapped.bam > PRJNA616446_mapped.bam
samtools sort PRJNA616446_mapped.bam
vi wuhCor1.genome
# NC_045512v2<TAB>29903
bedtools genomecov -d -ibam PRJNA616446_mapped.bam -g wuhCor1.genome > PRJNA616446_mapped.cov