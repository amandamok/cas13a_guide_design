##################################################
### compute folding energies along whole genome

library(optparse)
library(here)

option_list = list(make_option(c("-i", "--input"), type="character",
                               default=file.path(here(), "ref_data", "wuhCor1.RNAfold.txt"),
                               help="path to RNAfold output for wuhCor1 genome"),
                   make_option(c("-w", "--window"), type="integer", default=20,
                               help="window size", metavar="integer"),
                   make_option(c("-o", "--out"), type="character", default=".",
                               help="output directory", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat("\n")
cat("RNAfold score: whole genome\n")

# read in RNAfold output
RNAfold <- readLines(opt$input)
RNAfold_seq <- RNAfold[2]
RNAfold_structure <- strsplit(RNAfold[3], split=" ")[[1]][1]

# pull window starts
if(file.exists(file.path(opt$out, "windows_tags.txt"))) {
  cat("- pulling window positions from windows_tags.txt\n")
  window_starts <- read.table(file.path(opt$out, "windows_tags.txt"), header=T)$start
} else {
  window_starts <- seq.int(ncol(covid_alignment)-opt$window-3)
}

# split whole genome structure into windows
windows_structure <- sapply(window_starts$start,
                            function(x) {
                              substr(RNAfold_structure, start=x, stop=x+opt$window-1)
                            })

# write structure to output
RNAeval_input <- rep(NA, 2*length(windows_structure))
RNAeval_input[2*(seq.int(length(windows_structure)))-1] <- read.table(file.path(opt$out, "windows_tags.txt"),
                                                                      header=T, stringsAsFactors=F)$seq
RNAeval_input[2*(seq.int(length(windows_structure)))] <- windows_structure
writeLines(RNAeval_input,
           con=file.path(opt$out, "windows_genomeStructure.txt"))

# evaluate structure folding energy
if(!file.exists(file.path(opt$out, "windows_RNAfold_genome.txt"))) {
  cat("- evaluating folding energies\n")
  system(paste("/usr/bin/RNAeval -i", file.path(opt$out, "windows_genomeStructure.txt"),
         ">", file.path(opt$out, "windows_RNAfold_genome.txt")))
}
