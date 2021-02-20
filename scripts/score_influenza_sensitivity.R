##################################################
### count bowtie alignments

library(optparse)
library(here)
library(seqinr)
library(foreach)

option_list <- list(make_option(c("-t", "--type"), type="character", default=NULL,
                                help="virus type", metavar="integer"),
                    make_option(c("-m", "--mismatch"), type="character", default=1,
                                help="number of mismatches allowed", metavar="integer"),
                    make_option(c("-w", "--window"), type="integer", default=20,
                                help="length of spacer", metavar="integer"),
                    make_option(c("-o", "--out"), type="character", default=".",
                                help="output directory", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(!(opt$type %in% c("influenzaA", "influenzaB"))) {
  cat("\nERROR: incorrect virus type specified")
  q(save="no")
}

cat("\nCalculating sensitivity to influenza strains\n")

calculate_sensitivity <- function(alignment, consensus_seq, consensus_indices,
                                  window_start, mismatch=1, window_size=20) {
  ## return fraction of hits, allowing for #mismatch
  # alignment: matrix; output from as.matrix(seqinr::read.alignment())
  # consensus_seq: character vector; output from consensus(seqinr::read.alignment())
  # consensus_indices: integer vector; indices corresponding to A/T/C/G nucleotides in consensus
  # window_start: integer vector; start index of target window (relative to consensus)
  # mismatch: integer; number of allowed mismatches
  # window: integer; length of spacer window
  window_indices <- seq(consensus_indices[window_start],
                        consensus_indices[window_start + window_size - 1])
  strain_windows <- as.matrix(alignment)[, window_indices]
  num_mismatches <- sapply(seq(nrow(strain_windows)),
                           function(x) {
                             sum(strain_windows[x,] != consensus_seq[window_indices])
                           })
  num_hits <- sum(num_mismatches <= mismatch)
  return(num_hits / nrow(alignment))
}

# load alignments
alignments <- lapply(1:8,
                     function(x) {
                       fname <- file.path(here(), "ref_data", "influenza",
                                          paste0(opt$type, "_human_segment", x, ".aln"))
                       seqinr::read.alignment(fname, format="clustal")
                     })
consensus_seq <- lapply(alignments, consensus)
consensus_indices <- lapply(consensus_seq, function(x) { which(x != "-") })
alignment_matrices <- lapply(alignments, as.matrix)

# read in windows
windows <- read.table(file.path(opt$out, "windows.txt"),
                      header=T, stringsAsFactors=F)

# calculate sensitivity
windows$segment_num <- as.numeric(sub("segment", "", windows$segment))
windows$sensitivity <- sapply(seq(nrow(windows)),
                              function(x) {
                                segment_num <- windows$segment_num[x]
                                calculate_sensitivity(alignment_matrices[[segment_num]],
                                                      consensus_seq[[segment_num]],
                                                      consensus_indices[[segment_num]],
                                                      windows$start[x],
                                                      opt$mismatch,
                                                      opt$window)
                              })

# output computed sensitivity
write.table(windows[, c("segment", "start", "strand", "sensitivity")],
            file=file.path(opt$out, "score_sensitivity.txt"),
            quote=F, sep="\t", row.names=F)