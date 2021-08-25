# compile scores ----------------------------------------------------------

library(here)
source(file.path(here(), "scripts", "helper.R"))

# read in windows
windows <- read.table("windows.txt", header=T, stringsAsFactors=F, sep="\t")

# read in RNAfold crRNA score
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "has_crRNA_hairpin", "has_hairpin")
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "crRNA_spacer_basepairs", "spacer_basepairs")
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "gRNA_MFE", "MFE")

# read in alignments to other human coronaviruses
windows <- add_column(windows, "score_human_CoV_specificity.txt", "specificity", "specificity")

# read in alignments to human transcriptome
windows <- add_column(windows, "alignment_cts_GRCh38_latest_rna.txt", "match_against_hg38", "ct")

# read in alignments to cow transcriptome
windows <- add_column(windows, "alignment_cts_ARS-UCD1_rna.txt", "match_against_bosTau9", "ct")

# additional annotations --------------------------------------------------

# load fasta file
pop1_fasta <- readLines("human_pop1.fasta")[-1]
pop1_fasta <- paste(pop1_fasta, collapse="")
pop1_fasta <- strsplit(pop1_fasta, split="")[[1]]
pop1_fasta <- data.frame(pos = seq_along(pop1_fasta) + 98117293 - 1,
                         base = pop1_fasta)

# label exon/introns
pop1_fasta$exon <- ifelse(pop1_fasta$base %in% c("A", "T", "C", "G"),
                          "exon", "intron")
pop1_fasta$chunk <- 0
pop1_fasta$chunk[1] <- 1
for(x in 2:nrow(pop1_fasta)) {
  if(pop1_fasta$exon[x] != pop1_fasta$exon[x-1]) {
    pop1_fasta$chunk[x] <- pop1_fasta$chunk[x-1] + 1
  } else {
    pop1_fasta$chunk[x] <- pop1_fasta$chunk[x-1]
  }
}
exon_chunks <- unique(subset(pop1_fasta, exon=="exon")$chunk)
intron_chunks <- unique(subset(pop1_fasta, exon=="intron")$chunk)
exon_intron_labels <- c(paste0("exon", seq_along(exon_chunks)),
                        paste0("intron", seq_along(intron_chunks)))
names(exon_intron_labels) <- c(paste0("exon", exon_chunks),
                               paste0("intron", intron_chunks))
pop1_fasta$label <- exon_intron_labels[with(pop1_fasta, paste0(exon, chunk))]

# load common SNPs
pop1_snp <- read.table("human_pop1_snp151common.tsv",
                       col.names=c("chr", "start", "rsID", "snp_type", "alleles", "freq"),
                       stringsAsFactors=F)
pop1_snp$major_allele <- sapply(pop1_snp$alleles,
                                function(x) {
                                  strsplit(x, split=",")[[1]][1]
                                })
pop1_snp$minor_allele <- sapply(pop1_snp$alleles,
                                function(x) {
                                  strsplit(x, split=",")[[1]][2]
                                })
pop1_snp$major_freq <- sapply(pop1_snp$freq,
                              function(x) {
                                as.numeric(strsplit(x, split=",")[[1]][1])
                              })
pop1_snp$minor_freq <- sapply(pop1_snp$freq,
                              function(x) {
                                as.numeric(strsplit(x, split=",")[[1]][2])
                              })
pop1_snp$MAF <- sapply(seq(nrow(pop1_snp)),
                       function(x) {
                         min(pop1_snp$major_freq[x], pop1_snp$minor_freq[x])
                       })

# annotate SNP MAF
pop1_fasta$MAF <- 1
pop1_fasta$MAF[match(pop1_snp$start, pop1_fasta$pos)] <- pop1_snp$MAF
pop1_fasta$minor_freq <- 1
pop1_fasta$minor_freq[match(pop1_snp$start, pop1_fasta$pos)] <- pop1_snp$minor_freq

# add extra annotations ---------------------------------------------------

windows$label <- sapply(windows$start,
                        function(x) {
                          tmp_window <- pop1_fasta$label[x:(x+19)]
                          return(paste(unique(tmp_window), collapse="_"))
                        })
windows$exon_only <- grepl("^exon", windows$label) & !grepl("_", windows$label)

windows$num_SNPs <- sapply(windows$start,
                           function(x) {
                             tmp_window <- pop1_fasta$MAF[x:(x+19)]
                             return(sum(tmp_window != 1))
                           })

# write output ------------------------------------------------------------

write.table(windows, file="human_pop1_summary.txt", quote=F, col.names=T, sep="\t")
