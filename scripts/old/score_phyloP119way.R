##################################################
### score PhyloP

library(optparse)

option_list <- list(make_option(c("-i", "--input"), type="character", 
                               default="~/covid-19/ref_data/wuhCor1.phyloP119way.txt", 
                               help="multiz alignment file name", metavar="character"),
                   make_option(c("-w", "--window"), type="integer", default=20,
                               help="window size", metavar="integer"),
                   make_option(c("-o", "--out"), type="character", default=".", 
                               help="output directory", metavar="character")) 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("\n")
cat("PhyloP conservation score: 119way\n")

# read in PhyloP scores
phylop <- readLines(opt$input)
phylop <- phylop[grepl("^[[:digit:]]", phylop)]
phylop <- data.frame(matrix(as.numeric(unlist(strsplit(phylop, split="\t"))), ncol=2, byrow=T))
colnames(phylop) <- c("position", "PhyloP")

# compute average PhyloP score per window
if(file.exists(file.path(opt$out, "windows_tags.txt"))) {
  cat("- pulling window positions from windows_tags.txt\n")
  window_starts <- read.table(file.path(opt$out, "windows_tags.txt"), header=T)$start 
} else {
  window_starts <- seq.int(ncol(covid_alignment)-opt$window-3)
}
phylop_scores <- sapply(window_starts,
                        function(x) {
                          mean(phylop$PhyloP[x:(x+opt$window)])
                        })

# output PhyloP scores
write.table(data.frame(start=window_starts, phyloP=phylop_scores),
            file=file.path(opt$out, "score_phyloP_119way.txt"),
            quote=F, sep="\t", row.names=F)