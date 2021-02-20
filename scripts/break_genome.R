##################################################
### break genome into overlapping windows

library(optparse)
library(Biostrings, quietly=T)
library(here)

option_list <- list(make_option(c("-i", "--input"), type="character",
                                default=file.path(here(), "ref_data/NC_045512v2.fa"),
                                help="genome .fa file name", metavar="character"),
                    make_option(c("-w", "--window"), type="integer", default=20,
                                help="window size", metavar="integer"),
                    make_option(c("-e", "--enzyme"), type="character", default=NULL, # "Cas13a" or "Cas12"
                                help="Cas enzyme type", metavar="character"),
                    make_option(c("-g", "--genome_strand", type="character", default="+",
                                  help="strandedness of viral genome", metavar="character")),
                    make_option(c("-s", "--strand"), type="character", default="+",
                                help="strand to target", metavar="character"),
                    make_option(c("-o", "--out"), type="character", default=".",
                                help="output directory", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Breaking genome into windows\n")

if(!("enzyme" %in% names(opt))) {
  cat("\nNO ENZYME SPECIFIED\n")
  q()
}

# read in genome sequence, process into continuous string
# genome_seq <- readLines(opt$input)
# genome_seq <- genome_seq[!grepl(">", genome_seq)]
# genome_seq <- paste0(genome_seq, collapse="")
genome_seq <- readDNAStringSet(opt$input)
genome_seq <- as.character(genome_seq)
if(length(genome_seq) > 1) { names(genome_seq) <- paste0("segment", seq_along(genome_seq))}
# genome_seq <- gsub("-", "", genome_seq)

# break genome into windows
# windows <- data.frame(start=seq.int(nchar(genome_seq)-opt$window-3),
#                       target=sapply(seq.int(nchar(genome_seq)-opt$window-3), # all windows that have 4 nt following
#                                     function(x) {
#                                       substr(genome_seq, start=x, stop=x+opt$window-1)
#                                     }),
#                       stringsAsFactors=F)
windows <- lapply(seq_along(genome_seq),
                  function(x) {
                    segment_seq <- genome_seq[x]
                    num_windows <- nchar(segment_seq) - opt$window - 3
                    data.frame(segment = names(genome_seq)[x],
                               start=seq.int(num_windows),
                               target = sapply(seq.int(num_windows),
                                               function(x) {
                                                 substr(segment_seq,
                                                        start = x,
                                                        stop = x + opt$window - 1)
                                               }))
                  })
windows <- data.frame(do.call(rbind, windows), stringsAsFactors=F)
opposite_strand <- ifelse(opt$genome_strand == "+", "-", "+")
if(opt$enzyme == "Cas13a") {
  if(opt$strand == opt$genome_strand) {
    windows$spacer <- gsub("T", "U", as.character(reverseComplement(DNAStringSet(windows$target))))
    windows$strand <- opt$genome_strand
  } else {
    windows$spacer <- gsub("T", "U", windows$target)
    windows$target <- as.character(reverseComplement(DNAStringSet(windows$target)))
    windows$strand <- opposite_strand
  }
} else { # Cas12, add spacers that target minus strand
  windows_minusStrand <- windows
  windows_minusStrand$target <- gsub("T", "U", windows$spacer)
  windows_minusStrand$spacer <- gsub("T", "U", windows$target)
  windows_minusStrand$strand <- "-"
  windows <- rbind(windows, windows_minusStrand, stringsAsFactors=F)
}

# apply filters for tag/antitag and/or PAM
if(opt$enzyme == "Cas13a") {
  # windows$antitag <- sapply(windows$start,
  #                           function(x) {
  #                             substr(genome_seq, start=x+opt$window, stop=x+opt$window+3)
  #                           })
  windows$antitag <- sapply(seq(nrow(windows)),
                            function(x) {
                              substr(genome_seq[windows$segment[x]],
                                     start = windows$start[x] + opt$window,
                                     stop = windows$start[x] + opt$window + 3)
                            })
  windows$antitag <- gsub("T", "U", windows$antitag)
} else {
  if(opt$enzyme == "Cas12") {
    windows$PAM <- NA
    windows$PAM[windows$strand=="+"] <- sapply(windows$start[windows$strand=="+"],
                                               function(x) {
                                                 PAM <- substr(genome_seq, start=x+opt$window, stop=x+opt$window+3)
                                                 return(as.character(reverseComplement(DNAString(PAM))))
                                               })
    windows$PAM[windows$strand=="-"] <- sapply(windows$start[windows$strand=="-"],
                                               function(x) {
                                                 substr(genome_seq, start=x-4, stop=x-1)
                                               })
    windows <- subset(windows, nchar(PAM)>=4)
    windows <- subset(windows, grepl("^TTT", windows$PAM))
  }
}

# compute %GC and %A
windows$GC_content <- sapply(windows$spacer,
                             function(x) {
                               window_freq <- alphabetFrequency(RNAString(x))
                               return((window_freq["C"]+window_freq["G"])/opt$window)
                             })
windows$A_content <- sapply(windows$spacer,
                            function(x) {
                              return(alphabetFrequency(RNAString(x))["A"]/opt$window)
                            })

# compute distance to 3' end
# windows$dist_to_3primeEnd <- nchar(genome_seq)-windows$start

# make fasta file of filtered windows
# windows_fasta <- rep(windows$target, each=2)
# windows_fasta[2*seq.int(nrow(windows))-1] <- paste0(">", windows_fasta[2*seq.int(nrow(windows))-1])
windows_fasta <- DNAStringSet(windows$target)
names(windows_fasta) <- with(windows, paste(segment, start, sep="_"))

# write window sequences to output
writeLines(as.character(windows$target), con=file.path(opt$out, "targets.txt"))
writeLines(as.character(windows$spacer), con=file.path(opt$out, "spacers.txt"))
# writeLines(windows_fasta, con=file.path(opt$out, "targets.fa"))
writeXStringSet(windows_fasta, filepath=file.path(opt$out, "targets.fa"))
write.table(windows, file=file.path(opt$out, "windows.txt"),
            quote=F, sep="\t", row.names=F)
