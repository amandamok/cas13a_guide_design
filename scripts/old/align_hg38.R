##################################################
### count bowtie alignments

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

# align windows with bowtie
cts_fname <- file.path(opt$out, paste0("windows_", opt$genome, "_mapped_cts.txt"))
if(!file.exists(cts_fname)) {
  cat(paste("- aligning windows to", opt$genome, "(bowtie)"))
  system(paste("/mnt/ingolialab/linux-x86_64/bin/bowtie -a -v 3 -p 10 -S --un", 
               file.path(opt$out, paste0("windows_", opt$genome, "_unmapped.fa")), # fasta file of unmapped windows
               "-f", file.path("~/covid-19/ref_data", opt$genome),  # path to bowtie index
               file.path(opt$out, "windows.fa"), # fname of windows fasta file
               ">", file.path(opt$out, paste0("windows_", opt$genome, "_mapped.sam")), # sam alignment file of mapped windows
               "2>", file.path(opt$out, paste0("windows_", opt$genome, "_mapped.bowtiestats"))))
  system(paste("cut -f1", file.path(opt$out, paste0("windows_", opt$genome, "_mapped.sam")),
               "| grep -v @ | sort | uniq -c >", cts_fname))
}

# read in alignments
unmapped <- read.table(file.path(opt$out, paste0("windows_", opt$genome, "_unmapped.fa")), stringsAsFactors=F)$V1
unmapped <- unmapped[!grepl(">", unmapped)]
mapped <- readLines(cts_fname)
mapped <- data.frame(t(sapply(readLines(cts_fname),
                              function(x) {
                                tmp <- strsplit(x, split=" ")[[1]]
                                tmp <- tmp[!(tmp=="")]
                                return(tmp)
                              })), stringsAsFactors=F)
rownames(mapped) <- NULL
colnames(mapped) <- c("count", "seq")
mapped$count <- as.numeric(mapped$count)

# format table of alignment counts
windows <- read.table(file.path(opt$out, "windows_tags.txt"), header=T, stringsAsFactors=F)
duplicated_windows <- unique(windows[duplicated(windows$seq)])
duplicated_cts <- sapply(duplicated_windows,
                         function(x) {
                           sum(windows==x)
                         })
window_cts <- data.frame(seq=windows$seq,
                         start=windows$start,
                         ct=sapply(windows,
                                   function(x) {
                                     if(x %in% unmapped) {
                                       return(0)
                                     } else {
                                       if(x %in% duplicated_windows) {
                                         return(mapped$count[match(x, mapped$seq)]/duplicated_cts[x])
                                       } else {
                                         return(mapped$count[match(x, mapped$seq)])
                                       }
                                     }
                                   }))
cat("- number of alignments:\n")
print(data.frame("number of alignments"=0:5,
                 count=sapply(0:5, function(x) sum(window_cts$ct==x))))

# write alignment counts to output
write.table(window_cts,
            file=file.path(opt$out, paste0("alignment_cts_", opt$genome, ".txt")),
            quote=F, sep="\t", row.names=F)