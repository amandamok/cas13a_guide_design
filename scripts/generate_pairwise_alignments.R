##################################################
### align strains to reference

library(optparse)
library(Biostrings)
library(foreach)
library(here)

option_list <- list(make_option(c("-g", "--genome"), type="character",
                                default=file.path(here(), "ref_data/NC_045512v2.fa"),
                                help="genome .fa file name", metavar="character"),
                    make_option(c("-i", "--input"), type="character",
                                default=file.path(here(), "ref_data/gisaid_cov2020_sequences.fasta"),
                                help="SARS-CoV-2 .fasta sequences", metavar="character"),
                    make_option(c("-n", "--num_cores"), type="integer", default=1, 
                                help="number of cores to parallelize over", metavar="character"),
                    make_option(c("-c", "--chunk_size"), type="integer", default=10,
                                help="number of genomes to align at a time", metavar="character"),
                    make_option(c("-o", "--output"), type="character",
                                default=file.path(here(), "ref_data", "gisaid_cov2020_alignment.txt"),
                                help="output filepath", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("- aligning strains to reference\n")

wuhCor1 <- readDNAStringSet(opt$genome)
strains <- readDNAStringSet(opt$input)

# align strains to wuhCor1
num_cores <- opt$num_cores
chunk_size <- opt$chunk_size
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
write.table(alignment, file=opt$output,
            quote=F, col.names=F)