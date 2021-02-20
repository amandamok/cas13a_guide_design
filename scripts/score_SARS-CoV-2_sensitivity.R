##################################################
### calculating sensitivity to SARS-CoV-2 strains

library(optparse)
library(here)

option_list <- list(make_option(c("-g", "--genome"), type="character", 
                                default=file.path(here(), "ref_data/NC_045512v2.fa"), 
                                help="genome .fa file name", metavar="character"),
                    make_option(c("-i", "--input"), type="character",
                                default=file.path(here(), "ref_data/gisaid_cov2020_alignment.txt"),
                                help="SARS-CoV-2 .fasta sequences", metavar="character"),
                    make_option(c("-w", "--window"), type="integer", default=20,
                                help="window size", metavar="integer"),
                    make_option(c("-o", "--out"), type="character", default=".", 
                                help="output directory", metavar="character")) 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("\nCalculating sensitivity to SARS-CoV-2 strains\n")

# read in genome sequence, process into continuous string
genome_seq <- readLines(opt$genome)
genome_seq <- genome_seq[!grepl(">", genome_seq)]
genome_seq <- paste0(genome_seq, collapse="")
genome_seq <- strsplit(genome_seq, split="")[[1]]

# align SARS-CoV-2 genomes to wuhCor1
if(!file.exists(opt$input)) {
  system(paste("Rscript", file.path(here(), "scripts/generate_pairwise_alignments.R")))
} else {
  strains <- matrix(unlist(strsplit(readLines(opt$input), split=" ")), ncol=length(genome_seq)+1, byrow=T)
  rownames(strains) <- strains[,1]
  strains <- strains[,-1]
}

# evaluate mismatch
mismatch <- sapply(seq_along(genome_seq),
                   function(x) {
                     strains[,x] != genome_seq[x]
                   })
rownames(mismatch) <- rownames(strains)
colnames(mismatch) <- paste0("pos_", seq.int(ncol(mismatch)))

# score windows
if(file.exists(file.path(opt$out, "windows.txt"))) {
  cat("- pulling window positions from windows.txt\n")
  window_starts <- read.table(file.path(opt$out, "windows.txt"), header=T)$start 
} else {
  window_starts <- seq.int(ncol(length(genome_seq))-opt$window-3)
}
alignment_scores <- data.frame(t(sapply(window_starts,
                                        function(x) {
                                          tmp_window <- mismatch[, x:(x+opt$window-1)]
                                          tmp_num_mismatch <- rowSums(tmp_window)
                                          mismatch_0 <- sum(tmp_num_mismatch==0)
                                          mismatch_1 <- sum(tmp_num_mismatch==1)
                                          mismatch_2 <- sum(tmp_num_mismatch>=2)
                                          return(c(mismatch_0, mismatch_1, mismatch_2))
                                        })))
colnames(alignment_scores) <- paste0("gisaid_mismatch_", 0:2)

# output alignment scores
write.table(data.frame(start=window_starts,
                       strand=read.table(file.path(opt$out, "windows.txt"), header=T)$strand,
                       alignment_scores,
                       sensitivity_0=alignment_scores$gisaid_mismatch_0/nrow(strains),
                       sensitivity_01=(alignment_scores$gisaid_mismatch_0 + alignment_scores$gisaid_mismatch_1)/nrow(strains)),
            file=file.path(opt$out, "score_gisaid_SARS-CoV-2.txt"),
            quote=F, sep="\t", row.names=F)
