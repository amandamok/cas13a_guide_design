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
actb_snp <- read.table("human_actb_snp153common.tsv",
                       comment.char="", header=T, stringsAsFactors=F)
actb_snp$chromStart <- actb_snp$chromStart + 1 # accounting for UCSC coordinates

# annotate positions overlapping with variants in dbSNP 153
actb_fasta$dbSNP153 <- 0
for(x in seq(nrow(actb_snp))) {
  tmp_start <- actb_snp$chromStart[x]
  tmp_end <- actb_snp$chromEnd[x]
  which_pos <- match(seq(tmp_start, tmp_end), actb_fasta$pos)
  actb_fasta$dbSNP153[which_pos] <- actb_fasta$dbSNP153[which_pos] + 1
}

# add extra annotations ---------------------------------------------------

windows$label <- sapply(windows$start,
                        function(x) {
                          tmp_window <- actb_fasta$label[x:(x+19)]
                          return(paste(unique(tmp_window), collapse="_"))
                        })
windows$exon_only <- grepl("^exon", windows$label) & !grepl("_", windows$label)

windows$dbSNP153 <- sapply(windows$start,
                           function(x) {
                             tmp_window <- actb_fasta$dbSNP153[x:(x+19)]
                             return(sum(tmp_window != 0))
                           })

# write output ------------------------------------------------------------

windows$guide_id <- paste0("actb_", seq(nrow(windows)))
windows$bed_chrEnd <- 5563902 - windows$start + 1
windows$bed_chrStart <- windows$bed_chrEnd - 20 + 1

windows_bed <- data.frame(chr="chr7",
                          chrStart=windows$bed_chrStart,
                          chrEnd = windows$bed_chrEnd,
                          guide_id = windows$guide_id)

write.table(windows, file="human_actb_summary.txt",
            quote=F, row.names=F, col.names=T, sep="\t")
write.table(windows_bed, file="human_actb_spacers.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

# select guides -----------------------------------------------------------

# starting: n = 36728
# 1. crRNA: good secondary structure; avoid antitag G (n = 3483)
selected_windows <- subset(windows,
                           has_crRNA_hairpin &
                             crRNA_spacer_basepairs == 0 &
                             !grepl("^G", antitag))
# 2. sensitivity: avoid variant in dbSNP 153 common variants (n = 2830)
selected_windows <- subset(selected_windows,
                           dbSNP153 == 0)
# 3. specificty: avoid targeting other human coronaviruses or transcripts (n = 548)
selected_windows <- subset(selected_windows,
                           specificity == 1 &
                             match_against_hg38 == 0)
# 4. position: only target exonic regions (n = 9)
selected_windows <- subset(selected_windows, exon_only)

selected_windows_bed <- subset(windows_bed,
                               windows_bed$guide_id %in% selected_windows$guide_id)
write.table(selected_windows_bed, file="human_actb_selected.bed",
            quote=F, row.names=F, col.names=F, sep="\t")