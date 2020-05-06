##################################################
### count blat alignments

library(optparse)
option_list = list(make_option(c("-g", "--genome"), type="character", 
                               default=NULL, 
                               help="bowtie index prefix", metavar="character"),
                   make_option(c("-o", "--out"), type="character", default=".", 
                               help="output directory", metavar="character")) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat("\n")
if(is.null(opt$genome)) {
  cat("ERROR: no reference genome specified")
  q(save="no")
} else {
  cat(paste("aligning against off-target:", opt$genome, "(bowtie) \n"))
}


# align windows with BLAT
pslx_fname <- file.path(opt$out, paste0("windows_blat_", opt$genome, ".pslx"))
if(!file.exists(pslx_fname)) {
  cat(paste("- aligning windows to", opt$genome, "(BLAT)"))
  system(paste("/mnt/ingolialab/linux-x86_64/bin/blat", 
               file.path("~/covid-19/ref_data", paste0(opt$genome, ".fa")), 
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