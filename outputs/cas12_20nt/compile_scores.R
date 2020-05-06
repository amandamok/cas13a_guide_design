##################################################
### compile scores

source("~/covid-19/scripts/helper.R")

# read in windows
windows <- read.table("windows.txt", header=T, stringsAsFactors=F, sep="\t")

# read in alignment score to SARS-CoV-2 strains
windows <- add_column(windows, "score_gisaid_SARS-CoV-2.txt", "covid19_mismatch_0", "gisaid_mismatch_0")
windows <- add_column(windows, "score_gisaid_SARS-CoV-2.txt", "covid19_mismatch_1", "gisaid_mismatch_1")
windows <- add_column(windows, "score_gisaid_SARS-CoV-2.txt", "covid19_mismatch_2", "gisaid_mismatch_2")
windows <- add_column(windows, "score_gisaid_SARS-CoV-2.txt", "sensitivity_0", "sensitivity_0")
windows <- add_column(windows, "score_gisaid_SARS-CoV-2.txt", "sensitivity_01", "sensitivity_01")

# read in RNAfold crRNA score
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "has_crRNA_hairpin", "has_hairpin")
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "crRNA_spacer_basepairs", "spacer_basepairs")

# read in RNAfold target score
windows <- add_column(windows, "score_RNAfold_target.txt", "target_basepairing_prob", "target_basepairing_propensity")

# read in alignments to: other human coronaviruses
windows <- add_column(windows, "score_human_CoV_specificity.txt", "specificity", "specificity")

# read in alignments to: hg38
windows <- add_column(windows, "alignment_cts_hg38.txt", "match_against_hg38", "ct")

# read in alignments to: bosTau9
windows <- add_column(windows, "alignment_cts_bosTau9.txt", "match_against_bosTau9", "ct")

# read in alignments to: cross-reactive (co-occuring pathogens)
# windows <- add_column(windows, "alignment_cts_cross_reactive.txt", "match_against_pathogens", "ct")

# write output
write.table(windows, file="cas12_results_summary.txt", quote=F, col.names=T, sep="\t")