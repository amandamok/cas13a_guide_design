##################################################
### compile scores

source("~/cas13a_guide_design/scripts/helper.R")

# read inw windows
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
windows <- add_column(windows, "score_RNAfold_crRNAs.txt", "gRNA_MFE", "MFE")

# read in RNAfold target score
windows <- add_column(windows, "score_RNAfold_target.txt", "target_basepairing_prob", "target_basepairing_propensity")

# read in alignments to: other human coronaviruses
windows <- add_column(windows, "score_human_CoV_specificity.txt", "specificity", "specificity")

# read in alignments to: hg38
windows <- add_column(windows, "alignment_cts_GRCh38_latest_rna.txt", "match_against_hg38", "ct")

# read in alignments to: bosTau9
windows <- add_column(windows, "alignment_cts_ARS-UCD1_rna.txt", "match_against_bosTau9", "ct")

# read in alignments to: cross-reactive (co-occuring pathogens)
windows <- add_column(windows, "score_crossreactive_DNA.txt", "offtarget_rate", "DNA_offtarget_rate")
windows <- add_column(windows, "score_crossreactive_DNA.txt", "offtargets", "DNA_offtargets")
windows <- add_column(windows, "score_crossreactive_RNA.txt", "offtarget_rate", "RNA_offtarget_rate")
windows <- add_column(windows, "score_crossreactive_RNA.txt", "offtargets", "RNA_offtargets")

# write output
write.table(windows, file="cas13a_results_summary.txt", quote=F, col.names=T, sep="\t")
