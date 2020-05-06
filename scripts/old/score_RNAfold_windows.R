##################################################
### compute window folding energies

library(optparse)

option_list <- list(make_option(c("-o", "--out"), type="character", default=".", 
                               help="output directory", metavar="character")) 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("\nRNAfold score: windows\n")

# fold windows
if(!file.exists(file.path(opt$out, "windows_RNAfold.txt"))) {
  cat("- computing folding energies (windows)\n")
  system(paste0("/usr/bin/RNAfold -i ", file.path(opt$out, "windows.txt"), 
                " --noPS --outfile=", file.path("windows_RNAfold_window.txt")))
}
if(!file.exists(file.path(opt$out, "windows_reverseComplement_RNAfold.txt"))) {
  cat("- computing folding energies (windows_reverseComplement)\n")
  system(paste0("/usr/bin/RNAfold -i ", file.path(opt$out, "windows_reverseComplement.txt"),
                " --noPS --outfile=", file.path("windows_reverseComplement_RNAfold_window.txt")))
}

# read in RNAfold output: windows
RNAfold <- data.frame(matrix(readLines(file.path(opt$out, "windows_RNAfold_window.txt")), 
                             ncol=2, byrow=T), stringsAsFactors=F)
colnames(RNAfold) <- c("sequence", "structure")
RNAfold$target_MFE <- as.numeric(sub(" ", "", sub("\\)", "", sub(".*\\(", "", RNAfold$structure))))
RNAfold$target_structure <- sub(" .*", "", RNAfold$structure)

# read in RNAfold output: windows_reverseComplement
RNAfold_guide <- data.frame(matrix(readLines(file.path(opt$out, "windows_reverseComplement_RNAfold_window.txt")), 
                                   ncol=2, byrow=T), stringsAsFactors=F)
colnames(RNAfold_guide) <- c("sequence", "structure")
RNAfold$guide_MFE <- as.numeric(sub(" ", "", sub("\\)", "", sub(".*\\(", "", RNAfold_guide$structure))))
RNAfold$guide_structure <- sub(" .*", "", RNAfold_guide$structure)

# add start index
RNAfold$start <- read.table(file.path(opt$out, "windows_tags.txt"), header=T)$start

# output RNAfold scores
write.table(RNAfold, 
            file=file.path(opt$out, "score_RNAfold_windows.txt"),
            quote=F, sep="\t", row.names=F)
