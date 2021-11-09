##################################################
### calculate sensitivity with multisegment genome

library(optparse)
library(here)
library(seqinr)
library(foreach)

option_list <- list(make_option(c("-f", "--file"), type="character", default=NULL,
                                help="filepath header to alignment files", metavar="integer"),
                    make_option(c("-t", "--type"), type="character", default=NULL,
                                help="input file"),
                    make_option(c("-s", "--segments"), type="integer", default=8,
                                help="number of segments in input genome", metavar="integer"),
                    make_option(c("-m", "--mismatch"), type="character", default=1,
                                help="number of mismatches allowed", metavar="integer"),
                    make_option(c("-w", "--window"), type="integer", default=20,
                                help="length of spacer", metavar="integer"),
                    make_option(c("-o", "--out"), type="character", default=".",
                                help="output directory", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$type)) {
  cat("\nERROR: input file type not specified")
  q(save="no")
}

cat("\nCalculating sensitivity to sampled genomes\n")

calculate_sensitivity <- function(alignment, consensus_seq, consensus_indices,
                                  window_start, mismatch=1, window_size=20) {
  ## return fraction of hits, allowing for #mismatch
  # alignment: matrix; output from as.matrix(seqinr::read.alignment())
  # consensus_seq: character vector; output from consensus(seqinr::read.alignment())
  # consensus_indices: integer vector; indices corresponding to A/T/C/G nucleotides in consensus
  # window_start: integer vector; start index of target window (relative to consensus)
  # mismatch: integer; number of allowed mismatches
  # window: integer; length of spacer window
  consensus_window <- consensus_indices[seq(window_start,
                                            window_start + window_size - 1)]
  consensus_window <- consensus_seq[consensus_window]
  aln_window_start <- consensus_indices[window_start]
  strain_windows <- t(sapply(seq(nrow(alignment)),
                             function(x) {
                               tmp_stop <- aln_window_start + window_size - 1
                               tmp_window <- alignment[x, aln_window_start:tmp_stop]
                               num_nt <- sum(tmp_window != "-")
                               while(num_nt != window_size & tmp_stop < ncol(alignment)) {
                                 tmp_stop <- min(tmp_stop + (window_size - num_nt),
                                                 ncol(alignment))
                                 tmp_window <- alignment[x, aln_window_start:tmp_stop]
                                 num_nt <- sum(tmp_window != "-")
                               }
                               tmp_window <- tmp_window[tmp_window != "-"]
                               if(length(tmp_window) < window_size) {
                                 tmp_window <- c(tmp_window,
                                                 rep("-", window_size - length(tmp_window)))
                               }
                               return(tmp_window)
                             }))
  num_mismatches <- sapply(seq(nrow(strain_windows)),
                           function(x) {
                             sum(strain_windows[x,] != consensus_window)
                           })
  num_hits <- sum(num_mismatches <= mismatch)
  return(num_hits / nrow(alignment))
}

# load alignments
if(opt$type=="clustal") {
  alignments <- lapply(1:opt$segments,
                       function(x) {
                         fname <- paste0(opt$file, x, ".aln")
                         seqinr::read.alignment(fname, format="clustal")
                       })
} else {
  alignments <- lapply(1:opt$segments,
                       function(x) {
                         fname <- paste0(opt$file, x, ".fasta")
                         Biostrings::readDNAStringSet(fname)
                       })
}

alignment_matrices <- lapply(alignments, as.matrix)
consensus_seq <- lapply(alignment_matrices, seqinr::consensus)
consensus_indices <- lapply(consensus_seq, function(x) { which(x != "-") })

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