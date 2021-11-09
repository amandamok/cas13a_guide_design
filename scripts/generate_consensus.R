##################################################
### generate consensus sequence from CLUSTAL

library(optparse)
library(Biostrings, quietly=T)
library(seqinr, quietly=T)

option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL,
                                help="filepath to input alignment file", metavar="character"),
                    make_option(c("-t", "--type"), type="character", default=NULL,
                                help="input type (`clustal` or `fasta`)", metavar="character"),
                    make_option(c("-n", "--name"), type="character", default="consensus",
                                help="name of sequence", metavar="character"),
                    make_option(c("-o", "--output"), type="character", default="consensus.fasta",
                                help="filepath to output fasta file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$type)) {
  cat("\nERROR: input file type not specified")
  q(save="no")
}

if(!(opt$type %in% c("clustal", "fasta"))) {
  cat("\nERROR: input file type is not recognized")
  q(save="no")
}

if(opt$type == "clustal") {
  alignment <- seqinr::read.alignment(opt$input, format="clustal")
} else {
  alignment <- Biostrings::readDNAStringSet(opt$input)
  alignment <- as.matrix(alignment)
}
consensus <- seqinr::consensus(alignment, method="majority")
consensus <- gsub("-", "", paste(consensus, collapse=""))
consensus <- Biostrings::DNAStringSet(consensus)
names(consensus) <- opt$name
Biostrings::writeXStringSet(consensus, opt$output)
