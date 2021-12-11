##################################################
### score specificity against other human coronaviruses

library(optparse)
library(here)

option_list <- list(make_option(c("-g", "--genome"), type="character", default="human_CoV",
                                help="bowtie index prefix", metavar="character"),
                    make_option(c("-e", "--enzyme"), type="character", default=NULL, # "Cas13a" or "Cas12"
                                help="Cas enzyme type", metavar="character"),
                    make_option(c("-m", "--mismatch"), type="integer", default=1,
                                help="number of mismatches allowed", metavar="integer"),
                    make_option(c("-n", "--num_human_CoV"), type="integer", default=6,
                                help="number of other human coronaviruses", metavar="integer"),
                    make_option(c("-o", "--out"), type="character", default=".",
                                help="output directory", metavar="character"),
                    make_option(c("-b", "--bowtie"), type="character", default=NULL,
                                help="path to bowtie", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

bowtie_path <- system("which bowtie", intern=T)
if(is.null(opt$bowtie)) {
  opt$bowtie <- bowtie_path
} else {
  if(length(bowtie_path)==0) {
    cat("ERROR: need to supply path to bowtie")
    q(save="no")
  }
}

cat("\nCalculating specificity against other human coronaviruses\n")

# align windows with bowtie
cts_fname <- file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.sam"))
if(!file.exists(cts_fname)) {
  cat(paste("- aligning windows to", opt$genome, "\n"))
  system(paste(bowtie_path,
               ifelse(opt$enzyme=="Cas13a", "--norc", ""), # for Cas13a, do not align to reverse complement
               "-k 50", # report up to 50 alignments
               "-v", opt$mismatch, # up to opt$mismatch mismatches allowed
               "-S", # output as .sam alignment file
               "--un", file.path(opt$out, paste0("bowtie_", opt$genome, "_unmapped.fa")), # fasta file of unmapped windows
               "-f", file.path(here(), "ref_data", opt$genome),  # path to bowtie index
               file.path(opt$out, "targets.fa"), # fname of targets fasta file
               ">", file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.sam")), # sam alignment file of mapped windows
               "2>", file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.bowtiestats")))) # fname of bowtie output
}

# read in alignment file
alignment <- system(paste("grep -v ^@", file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.sam")),
                          "| cut -f1,3"), intern=T)
if(length(alignment)==0) {
  cat("- no alignments reported\n")
} else {
  alignment <- data.frame(matrix(unlist(strsplit(alignment, split="\t")), ncol=2, byrow=T), stringsAsFactors=F)
  colnames(alignment) <- c("window", "aligned_to")
  alignment <- subset(alignment, alignment$aligned_to != "*")
}

# compute specificity
windows <- read.table(file.path(opt$out, "windows.txt"), header=T, stringsAsFactors=F)
if(length(alignment)==0) {
  num_CoV_aligned <- 0
} else {
  num_CoV_aligned <- sapply(windows$target,
                            function(x) {
                              if(!(x %in% alignment$window)) {
                                return(0)
                              } else {
                                tmp <- subset(alignment, alignment$window==x)
                                return(length(unique(tmp$aligned_to)))
                              }
                            })
}
# num_human_CoV <- as.numeric(system(paste("grep ^'>'", file.path(here(), "ref_data", paste0(opt$genome, ".fa")),
#                                          "| wc -l"), intern=T))
specificity <- (opt$num_human_CoV - num_CoV_aligned) / opt$num_human_CoV

# output specificity score
write.table(data.frame(segment=windows$segment,
                       start=windows$start,
                       strand=windows$strand,
                       num_CoV_aligned, specificity),
            file=file.path(opt$out, "score_human_CoV_specificity.txt"),
            quote=F, sep="\t", row.names=F)
