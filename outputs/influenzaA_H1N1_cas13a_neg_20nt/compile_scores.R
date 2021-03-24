##################################################
### compile scores

source("../../scripts/helper.R")

# read in windows
windows <- read.table("windows.txt", header=T, stringsAsFactors=F, sep="\t")

# read in sensitivity scores
windows <- add_column(windows, "score_sensitivity.txt", "sensitivity", "sensitivity")

# read in RNAfold crRNA score
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "has_crRNA_hairpin", "has_hairpin")
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "crRNA_spacer_basepairs", "spacer_basepairs")

# read in alignments to: hg38
windows <- add_column(windows, "alignment_cts_GRCh38_latest_rna.txt", "match_against_hg38", "ct")

# read in alignments to: bosTau9
windows <- add_column(windows, "alignment_cts_ARS-UCD1_rna.txt", "match_against_bosTau9", "ct")

# write output
write.table(windows, file="influenzaA_H1N1_cas13a_neg_results_summary.txt", quote=F, col.names=T, sep="\t")
