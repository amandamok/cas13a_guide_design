##################################################
### count bowtie alignments

library(optparse)
library(ggplot2)

option_list <- list(make_option(c("-g", "--genome"), type="character", 
                                default="~/covid-19/ref_data/wuhCor1", 
                                help="genome index prefix", metavar="character"),
                    make_option(c("-a", "--accession"), type="character",
                                default=NULL,
                                help="accession number to .fastq file", metavar="character"),
                    make_option(c("-p", "--paired"), action="store_true", default=F,
                                help="data is paired-end"),
                    make_option(c("-o", "--out"), type="character",
                                default="~/covid-19/ref_data/RNA_expression",
                                help="output directory", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(!("accession" %in% names(opt))) {
  cat("\nNO DATASET SPECIFIED\n")
  q()
}

cat("\n- aligning to wuhCor1 genome\n")

alignment_fname <- file.path(opt$out, paste0(opt$accession, "_mapped.sam"))
if(!file.exists(alignment_fname)) {
  if(!file.exists(file.path(opt$out, paste0(opt$accession, ".fastq")))) {
    system(paste("prefetch", opt$accession))
    if(opt$paired) {
      system(paste("fastq-dump --split-files -O", opt$out, file.path("~/ncbi/public/sra", paste0(opt$accession, ".sra"))))
      system(paste("cat", file.path(opt$out, paste0(opt$accession, "_1.fastq")), 
                   file.path(opt$out, paste0(opt$accession, "_2.fastq")),
                   ">", file.path(opt$out, paste0(opt$accession, ".fastq"))))
    } else {
      system(paste("fastq-dump", "-O", opt$out, file.path("~/ncbi/public/sra", paste0(opt$accession, ".sra"))))
    }
  }
  system(paste("/mnt/ingolialab/linux-x86_64/bin/bowtie -S -q ~/covid-19/ref_data/wuhCor1", 
               file.path(opt$out, paste0(opt$accession, ".fastq")),
               ">", file.path(opt$out, paste0(opt$accession, "_mapped.sam")),
               "2>", file.path(opt$out, paste0(opt$accession, "_mapped.bowtiestats"))))
}

aligned_reads <- system(paste("grep NC_045512v2", file.path(opt$out, paste0(opt$accession, "_mapped.sam")),
                              "| grep -v ^@ | cut -f1,3,4"), intern=T)
aligned_reads <- data.frame(matrix(unlist(strsplit(aligned_reads, split="\t")), ncol=3, byrow=T), stringsAsFactors=F)
colnames(aligned_reads) <- c("seqID", "rname", "pos")
aligned_reads$pos <- as.numeric(aligned_reads$pos)

pileup_plot <- ggplot(aligned_reads, aes(pos)) + geom_histogram(binwidth=300) + 
  theme_bw() + ggtitle(opt$accession) + xlab("genomic position") + ylab("read count")
ggsave(plot=pileup_plot,
       filename=file.path(opt$out, paste0(opt$accession, "_plot.pdf")),
       device="pdf", width=8, height=6, units="in")

