##################################################
### count blat alignments

library(optparse)
library(here)

option_list = list(make_option(c("-g", "--genome"), type="character",
                               default=NULL,
                               help="bowtie index prefix", metavar="character"),
                   make_option(c("-o", "--out"), type="character", default=".",
                               help="output directory", metavar="character"),
                   make_option(c("-b", "--blat"), type="character", default=NULL,
                               help="path to blat", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

blat_path <- system("which blat", intern=T)
if(length(blat_path)==0 & is.null(opt$blat)) {
  cat("ERROR: need to supply path to blat")
  q(save="no")
}

cat("\n")
if(is.null(opt$genome)) {
  cat("ERROR: no reference genome specified")
  q(save="no")
}

# align windows with BLAT
pslx_fname <- file.path(opt$out, paste0("windows_blat_", opt$genome, ".pslx"))
if(!file.exists(pslx_fname)) {
  cat(paste("- aligning windows to", opt$genome, "(BLAT)"))
  system(paste(blat_path,
               file.path(here(), "ref_data", paste0(opt$genome, ".fa")),
               file.path(opt$out, "windows.fa"),
               "-out=pslx", pslx_fname))
}

# read in BLAT output
blat <- readLines(pslx_fname)
blat <- blat[-c(1,2)] # remove header lines
blat <- blat[-length(blat)] # remove footer line
blat <- strsplit(blat, split="\t")
if(length(blat)==2) {
  cat("- no BLAT alignments")
  q(save="no")
}