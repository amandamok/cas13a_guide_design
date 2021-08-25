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
actb_fasta <- readLines("human_actb.fasta")[-1]
actb_fasta <- paste(actb_fasta, collapse="")
actb_fasta <- strsplit(actb_fasta, split="")[[1]]
actb_fasta <- data.frame(pos = rev(seq_along(actb_fasta) + 5527152 - 1),
                         base = actb_fasta,
                         stringsAsFactors=F)

# label exon/introns
actb_fasta$exon <- ifelse(actb_fasta$base %in% c("A", "T", "C", "G"),
                          "exon", "intron")
actb_fasta$chunk <- 0
actb_fasta$chunk[1] <- 1
for(x in 2:nrow(actb_fasta)) {
  if(actb_fasta$exon[x] != actb_fasta$exon[x-1]) {
    actb_fasta$chunk[x] <- actb_fasta$chunk[x-1] + 1
  } else {
    actb_fasta$chunk[x] <- actb_fasta$chunk[x-1]
  }
}
exon_chunks <- unique(subset(actb_fasta, exon=="exon")$chunk)
intron_chunks <- unique(subset(actb_fasta, exon=="intron")$chunk)
exon_intron_labels <- c(paste0("exon", seq_along(exon_chunks)),
                        paste0("intron", seq_along(intron_chunks)))
names(exon_intron_labels) <- c(paste0("exon", exon_chunks),
                               paste0("intron", intron_chunks))
actb_fasta$label <- exon_intron_labels[with(actb_fasta, paste0(exon, chunk))]

# load common SNPs
actb_snp <- read.table("human_actb_snp151common.tsv",
                       col.names=c("chr", "start", "rsID", "snp_type", "alleles", "freq"),
                       stringsAsFactors=F)
actb_snp$major_allele <- sapply(actb_snp$alleles,
                                function(x) {
                                  strsplit(x, split=",")[[1]][1]
                                })
actb_snp$minor_allele <- sapply(actb_snp$alleles,
                                function(x) {
                                  strsplit(x, split=",")[[1]][2]
                                })
actb_snp$major_freq <- sapply(actb_snp$freq,
                              function(x) {
                                as.numeric(strsplit(x, split=",")[[1]][1])
                              })
actb_snp$minor_freq <- sapply(actb_snp$freq,
                              function(x) {
                                as.numeric(strsplit(x, split=",")[[1]][2])
                              })
actb_snp$MAF <- sapply(seq(nrow(actb_snp)),
                       function(x) {
                         min(actb_snp$major_freq[x], actb_snp$minor_freq[x])
                       })

# annotate SNP MAF
actb_fasta$MAF <- 1
actb_fasta$MAF[match(actb_snp$start, actb_fasta$pos)] <- actb_snp$MAF
actb_fasta$minor_freq <- 1
actb_fasta$minor_freq[match(actb_snp$start, actb_fasta$pos)] <- actb_snp$minor_freq

# add extra annotations ---------------------------------------------------

windows$label <- sapply(windows$start,
                        function(x) {
                          tmp_window <- actb_fasta$label[x:(x+19)]
                          return(paste(unique(tmp_window), collapse="_"))
                        })
windows$exon_only <- grepl("^exon", windows$label) & !grepl("_", windows$label)

windows$num_SNPs <- sapply(windows$start,
                           function(x) {
                             tmp_window <- actb_fasta$MAF[x:(x+19)]
                             return(sum(tmp_window != 1))
                           })

# write output ------------------------------------------------------------

write.table(windows, file="human_actb_summary.txt", quote=F, col.names=T, sep="\t")
