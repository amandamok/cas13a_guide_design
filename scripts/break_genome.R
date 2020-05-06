##################################################
### break genome into overlapping windows

library(optparse)
library(Biostrings, quietly=T)

option_list <- list(make_option(c("-i", "--input"), type="character", 
                                default="~/covid-19/ref_data/NC_045512v2.fa", 
                                help="genome .fa file name", metavar="character"),
                    make_option(c("-w", "--window"), type="integer", default=20,
                                help="window size", metavar="integer"),
                    make_option(c("-e", "--enzyme"), type="character", default=NULL, # "Cas13a" or "Cas12"
                                help="Cas enzyme type", metavar="character"),
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
genome_seq <- readLines(opt$input)
genome_seq <- genome_seq[!grepl(">", genome_seq)]
genome_seq <- paste0(genome_seq, collapse="")

# break genome into windows
windows <- data.frame(start=seq.int(nchar(genome_seq)-opt$window-3),
                      target=sapply(seq.int(nchar(genome_seq)-opt$window-3), # all windows that have 4 nt following
                                    function(x) {
                                      substr(genome_seq, start=x, stop=x+opt$window-1)
                                    }),
                      stringsAsFactors=F)
windows$spacer <- gsub("T", "U", as.character(reverseComplement(DNAStringSet(windows$target))))
windows$strand <- "+"
if(opt$enzyme == "Cas12") { # also target minus strand of dsDNA
  windows_minusStrand <- windows
  windows_minusStrand$target <- gsub("T", "U", windows$spacer)
  windows_minusStrand$spacer <- gsub("T", "U", windows$target)
  windows_minusStrand$strand <- "-"
  windows <- rbind(windows, windows_minusStrand, stringsAsFactors=F)
}

# apply filters for tag/antitag and/or PAM
if(opt$enzyme == "Cas13a") {
  windows$antitag <- next_window <- sapply(windows$start,
                                           function(x) {
                                             substr(genome_seq, start=x+opt$window, stop=x+opt$window+3)
                                           })
  windows$antitag <- gsub("T", "U", windows$antitag)
  # windows <- subset(windows, !(windows$antitag=="GTTT"))
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
windows$dist_to_3primeEnd <- nchar(genome_seq)-windows$start

# make fasta file of filtered windows
windows_fasta <- rep(windows$target, each=2)
windows_fasta[2*seq.int(nrow(windows))-1] <- paste0(">", windows_fasta[2*seq.int(nrow(windows))-1])

# write window sequences to output
writeLines(windows$target, con=file.path(opt$out, "targets.txt"))
writeLines(windows$spacer, con=file.path(opt$out, "spacers.txt"))
writeLines(windows_fasta, con=file.path(opt$out, "targets.fa"))
write.table(windows, file=file.path(opt$out, "windows.txt"),
            quote=F, sep="\t", row.names=F)
