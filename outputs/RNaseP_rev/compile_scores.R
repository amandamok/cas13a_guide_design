##################################################
### compile scores

source("~/covid-19/scripts/helper.R")

# read in windows
windows <- read.table("windows.txt", header=T, stringsAsFactors=F, sep="\t")

# read in RNAfold crRNA score
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "has_crRNA_hairpin", "has_hairpin")
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "crRNA_spacer_basepairs", "spacer_basepairs")

# read in alignments to: hg38
windows <- add_column(windows, "alignment_cts_GRCh38_latest_rna.txt", "match_against_hg38", "ct")

# read in alignments to: wuhCor1
windows <- add_column(windows, "alignment_cts_wuhCor1.txt", "match_against_wuhCor1", "ct")

# write output
write.table(windows, file="RNaseP_rev_summary.txt", quote=F, col.names=T, sep="\t")
