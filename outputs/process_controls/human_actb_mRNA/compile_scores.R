library(here)
library(Biostrings)

source(file.path(here(), "scripts", "helper.R"))

# https://www.ncbi.nlm.nih.gov/nuccore/NG_007992.1?from=4482&to=8971&report=graph&content=5

# generate scores ---------------------------------------------------------

# concatenate exons -> generate fasta
exons_seq <- readDNAStringSet("NM_001101.3.exons.fasta")
mature_mRNA <- paste0(as.character(exons_seq), collapse="")
writeLines(c(">NM_001101.3", mature_mRNA),
           "mature_mRNA.fasta")

# generate scores
system("./process.sh")

# compile scores ----------------------------------------------------------

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

# load gff
gff3 <- readLines("Genes.GFF3")
gff3 <- gff3[!grepl("^#", gff3)]
gff3 <- data.frame(matrix(unlist(strsplit(gff3, split="\t")),
                          ncol=9, byrow=T),
                   stringsAsFactors=F)
colnames(gff3) <- c("seqid", "source", "type", "start", "end", "score",
                    "strand", "phase", "attributes")
gff3$start <- as.numeric(gff3$start)
gff3$end <- as.numeric(gff3$end)
gff3$width <- with(gff3, end - start + 1)
exons <- subset(gff3, grepl("exon.*number", gff3$attributes))

# link mature mRNA position to genomic position
mature_mRNA <- data.frame(position=seq(nchar(mature_mRNA)),
                          nt=strsplit(mature_mRNA, split="")[[1]])
mature_mRNA$index <- NA
mature_mRNA$exon_label <- NA
for(x in seq(nrow(exons))) {
  if(x == 1) {
    pos_start <- 1
  } else {
    pos_start <- sum(exons$width[1:(x-1)]) + 1
  }
  pos_end <- pos_start + exons$width[x] - 1
  mature_mRNA$index[seq(pos_start, pos_end)] <- seq(exons$start[x], exons$end[x])
  mature_mRNA$exon_label[seq(pos_start, pos_end)] <- paste0("exon", x)
}

# annotate SNPs
dbSNP_155 <- read.table("dbSNP_b155_v2.BED", stringsAsFactors=F)
colnames(dbSNP_155) <- c("name", "start", "end", "SNP", "score", "phase")
dbSNP_155 <- subset(dbSNP_155, end %in% mature_mRNA$index)
mature_mRNA$dbSNP_155 <- NA
mature_mRNA$dbSNP_155[match(dbSNP_155$end, mature_mRNA$index)] <- dbSNP_155$SNP

# add extra annotations ---------------------------------------------------

windows$label <- sapply(windows$start,
                        function(x) {
                          tmp_window <- mature_mRNA$exon_label[x:(x+19)]
                          return(paste(unique(tmp_window), collapse="_"))
                        })

windows$dbSNP155 <- sapply(windows$start,
                           function(x) {
                             tmp_window <- mature_mRNA$dbSNP_155[x:(x+19)]
                             return(sum(!is.na(tmp_window)))
                           })

# write output ------------------------------------------------------------

windows$guide_id <- paste0("actb_mRNA_", seq(nrow(windows)))

write.table(windows, file="human_actb_mRNA_summary.txt",
            quote=F, row.names=F, col.names=T, sep="\t")

# select guides -----------------------------------------------------------

# starting: n = 1789
# 1. crRNA: good secondary structure; avoid antitag G (n = 182)
selected_windows <- subset(windows,
                           has_crRNA_hairpin &
                             crRNA_spacer_basepairs == 0 &
                             !grepl("^G", antitag))
# 2. sensitivity: avoid variant in dbSNP 153 common variants (n = 174)
selected_windows <- subset(selected_windows,
                           dbSNP155 == 0)
# 3. specificty: avoid targeting other human coronaviruses or transcripts (n = 11)
selected_windows <- subset(selected_windows,
                           specificity == 1 &
                             match_against_hg38 == 0)

previous_set <- read.csv("../new_process_controls.csv")
new_set <- subset(selected_windows, !(target %in% previous_set$target))
## new set targets 2 SNPs in dbSNP_153: rs7612 and rs11546906