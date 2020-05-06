##################################################
### compile scores

library(optparse)
library(prodlim)

option_list = list(make_option(c("-o", "--out"), type="character", default=".", 
                               help="output directory", metavar="character")) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat("\nCompiling scores\n")

add_column <- function(dat, fname, dat_column, fname_column) {
  # add score from fname to dat
  ## dat: data.frame; contains columns called "start" and "strand
  ## fname: character; file.path to score to be added, file should include columns called "start" and "strand
  ## dat_column: character; name of column in dat to be added/replaced
  ## fname_column: character; name of column in fname to be pulled
  if(!file.exists(fname)) {
    stop(paste(fname, "does not exist"))
  }
  score_dat <- read.table(fname, header=T, stringsAsFactors=F, sep="\t")
  row_indices <- row.match(dat[, c("start", "strand")], score_dat[, c("start", "strand")])
  dat[dat_column] <- score_dat[row_indices, fname_column]
  return(dat)
}

# read in windows and tags
windows <- read.table("windows.txt", header=T, stringsAsFactors=F, sep="\t")

# read in alignment score to SARS-CoV-2 strains
windows <- add_column(windows, file.path(opt$out, "score_gisaid_SARS-CoV-2.txt"),
                      "covid19_mismatch_0", "gisaid_mismatch_0")
windows <- add_column(windows, file.path(opt$out, "score_gisaid_SARS-CoV-2.txt"),
                      "covid19_mismatch_1", "gisaid_mismatch_1")
windows <- add_column(windows, file.path(opt$out, "score_gisaid_SARS-CoV-2.txt"),
                      "covid19_mismatch_2", "gisaid_mismatch_2")
windows <- add_column(windows, file.path(opt$out, "score_gisaid_SARS-CoV-2.txt"),
                      "sensitivity_0", "sensitivity_0")
windows <- add_column(windows, file.path(opt$out, "score_gisaid_SARS-CoV-2.txt"),
                      "sensitivity_01", "sensitivity_01")

# read in RNAfold crRNA score
windows <- add_column(windows, file.path(opt$out, "score_RNAfold_crRNAs.txt"),
                      "has_crRNA_hairpin", "has_hairpin")
windows <- add_column(windows, file.path(opt$out, "score_RNAfold_crRNAs.txt"),
                      "crRNA_spacer_basepairs", "spacer_basepairs")

# read in RNAfold target score
windows <- add_column(windows, file.path(opt$out, "score_RNAfold_target.txt"),
                      "target_basepairing_prob", "target_basepairing_propensity")

# read in alignments to: other human coronaviruses
windows <- add_column(windows, file.path(opt$out, "score_human_CoV_specificity.txt"),
                      "specificity", "specificity")

# read in alignments to: hg38
windows <- add_column(windows, file.path(opt$out, "alignment_cts_hg38.txt"),
                      "match_against_hg38", "ct")

# read in alignments to: bosTau9
windows <- add_column(windows, file.path(opt$out, "alignment_cts_bosTau9.txt"),
                      "match_against_bosTau9", "ct")

# read in alignments to: cross-reactive (co-occuring pathogens)
windows <- add_column(windows, file.path(opt$out, "alignment_cts_cross_reactive.txt"),
                      "match_against_pathogens", "ct")

write.table(windows, file=file.path(opt$out, "results_summary.txt"),
            quote=F, col.names=T, sep="\t")