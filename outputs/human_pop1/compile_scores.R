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
pop1_snp <- read.table("human_pop1_snp153common.tsv", sep="\t",
                       comment.char="", header=T, stringsAsFactors=F)
pop1_snp$chromStart <- pop1_snp$chromStart + 1 # accounting for UCSC coordinates

# annotate positions overlapping with variants in dbSNP 153
pop1_fasta$dbSNP153 <- 0
for(x in seq(nrow(pop1_snp))) {
  tmp_start <- pop1_snp$chromStart[x]
  tmp_end <- pop1_snp$chromEnd[x]
  which_pos <- match(seq(tmp_start, tmp_end), pop1_fasta$pos)
  pop1_fasta$dbSNP153[which_pos] <- pop1_fasta$dbSNP153[which_pos] + 1
}

# add extra annotations ---------------------------------------------------

windows$label <- sapply(windows$start,
                        function(x) {
                          tmp_window <- pop1_fasta$label[x:(x+19)]
                          return(paste(unique(tmp_window), collapse="_"))
                        })
windows$exon_only <- grepl("^exon", windows$label) & !grepl("_", windows$label)

windows$dbSNP153 <- sapply(windows$start,
                           function(x) {
                             tmp_window <- pop1_fasta$dbSNP153[x:(x+19)]
                             return(sum(tmp_window != 0))
                           })

# write output ------------------------------------------------------------

windows$guide_id <- paste0("pop1_", seq(nrow(windows)))
windows$bed_chrStart <- 98117293 + windows$start - 1
windows$bed_chrEnd <- windows$bed_chrStart + 20 - 1

windows_bed <- data.frame(chr="chr7",
                          chrStart=windows$bed_chrStart,
                          chrEnd = windows$bed_chrEnd,
                          guide_id = windows$guide_id)

write.table(windows, file="human_pop1_summary.txt",
            quote=F, row.names=F, col.names=T, sep="\t")
write.table(windows_bed, file="human_pop1_spacers.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

# select guides -----------------------------------------------------------

# starting: n = 42520
# 1. crRNA: good secondary structure; avoid antitag G (n = 9843)
selected_windows <- subset(windows,
                           has_crRNA_hairpin &
                             crRNA_spacer_basepairs == 0 &
                             !grepl("^G", antitag))
# 2. sensitivity: avoid variant in dbSNP 153 common variants (n = 8589)
selected_windows <- subset(selected_windows,
                           dbSNP153 == 0)
# 3. specificty: avoid targeting other human coronaviruses or transcripts (n = 3351)
selected_windows <- subset(selected_windows,
                           specificity == 1 &
                             match_against_hg38 == 0)

selected_windows_bed <- subset(windows_bed,
                               windows_bed$guide_id %in% selected_windows$guide_id)
write.table(selected_windows_bed, file="human_pop1_selected.bed",
            quote=F, row.names=F, col.names=F, sep="\t")
