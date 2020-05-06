##################################################
### align SARS-CoV-2 strains to wuhCor1 reference

library(optparse)
library(Biostrings)
library(foreach)

option_list <- list(make_option(c("-g", "--genome"), type="character", 
                                default="~/covid-19/ref_data/NC_045512v2.fa", 
                                help="genome .fa file name", metavar="character"),
                    make_option(c("-i", "--input"), type="character",
                                default="~/covid-19/ref_data/gisaid_cov2020_sequences.fasta",
                                help="SARS-CoV-2 .fasta sequences", metavar="character"),
                    make_option(c("-o", "--out"), type="character",
                                default="~/covid-19/ref_data",
                                help="output directory", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("- aligning SARS-CoV-2 strains to wuhCor1 reference\n")

wuhCor1 <- readDNAStringSet(opt$genome)
strains <- readDNAStringSet(opt$input)

# align strains to wuhCor1
num_cores <- 15
chunk_size <- 10
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)
strains_chunk <- split(strains, ceiling(seq_along(strains)/chunk_size))
alignment <- foreach(x=seq_along(strains_chunk), 
                     .combine='rbind', .packages="Biostrings", .inorder=F) %dopar% {
                       tmp_align <- pairwiseAlignment(strains_chunk[[x]], wuhCor1)
                       tmp_align <- data.frame(as.matrix(tmp_align), stringsAsFactors=F)
                       rownames(tmp_align) <- names(strains_chunk[[x]])
                       return(tmp_align)
                     }
parallel::stopCluster(cl)

# write alignment output
output_fname <- "gisaid_cov2020_alignment.txt"
write.table(alignment, file=file.path(opt$out, output_fname),
            quote=F, col.names=F)