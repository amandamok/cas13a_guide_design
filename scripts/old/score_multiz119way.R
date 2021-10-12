##################################################
### score multiz alignment

library(optparse)
library(here)

option_list = list(make_option(c("-i", "--input"), type="character",
                               default=file.path(here(), "ref_data", "wuhCor1.multiz119way.txt"),
                               help="multiz alignment file name", metavar="character"),
                   make_option(c("-w", "--window"), type="integer", default=20,
                               help="window size", metavar="integer"),
                   make_option(c("-o", "--out"), type="character", default=".",
                               help="output directory", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat("\n")
cat("Multiz alignment score: among SARS-CoV-2 strains\n")

# read in multiz alignment
multiz <- readLines(opt$input)
multiz <- multiz[grepl("^s", multiz)]
multiz <- data.frame(t(sapply(multiz,
                              function(x) {
                                tmp <- strsplit(x, split=" ")[[1]]
                                tmp <- tmp[tmp != ""]
                                names(tmp) <- NULL
                                return(tmp[2:7])
                              })), stringsAsFactors=F)
rownames(multiz) <- NULL
colnames(multiz) <- c("src", "start", "size", "strand", "srcSize", "text")

# pull sequence of only COVID-19 sequences
covid <- subset(multiz, !grepl("^NC_", multiz$src))
covid_alignment <- t(sapply(covid$text,
                            function(x) {
                              strsplit(x, split="")[[1]]
                            }))
rownames(covid_alignment) <- covid$src

# remove positions where wuhCor1 alignment sequence is "-"
which_wuhCor1 <- grep("wuhCor1", rownames(covid_alignment))
covid_gaps <- subset(covid_alignment,
                     select=(covid_alignment[which_wuhCor1,]=="-"))
cat("- number of positions where gap in wuhCor1 is not gap in other SARS-CoV-2 alignments\n")
summary(as.factor(colSums(covid_gaps!="-")))
non_gapped <- sapply(seq.int(ncol(covid_gaps)),
                     function(x) {
                       if(any(covid_gaps[,x]!="-")) {
                         return(rownames(covid_gaps)[covid_gaps[,x]!="-"])
                       } else {
                         return(NA)
                       }
                     })
cat("- SARS-CoV-2 alignments where wuhCor1 is gapped\n")
summary(as.factor(non_gapped))
covid_alignment <- subset(covid_alignment,
                          select=(covid_alignment[which_wuhCor1,]!="-"))

# check that wuhCor1 alignment sequence corresponds to genome fasta file
genome_seq <- readLines(file.path(here(), "ref_data", "NC_045512v2.fa"))
genome_seq <- genome_seq[!grepl(">", genome_seq)]
genome_seq <- paste0(genome_seq, collapse="")
wuhCor1_seq <- paste0(covid_alignment[which_wuhCor1,], collapse="")
cat("- wuhCor1 alignment sequence corresponds to genome .fasta sequence\n")
cat(paste(" ", genome_seq == wuhCor1_seq, "\n"))

# score windows
if(file.exists(file.path(opt$out, "windows_tags.txt"))) {
  cat("- pulling window positions from windows_tags.txt\n")
  window_starts <- read.table(file.path(opt$out, "windows_tags.txt"), header=T)$start
} else {
  window_starts <- seq.int(ncol(covid_alignment)-opt$window-3)
}
alignment_scores <- sapply(window_starts,
                           function(x) {
                             pos_scores <- sapply(x:(x+opt$window),
                                                  function(y) {
                                                    mean(covid_alignment[-which_wuhCor1, y] ==
                                                           covid_alignment[which_wuhCor1, y])
                                                  })
                             mean(pos_scores)
                           })

# output alignment scores
write.table(data.frame(start=window_starts, multiz_score=alignment_scores),
            file=file.path(opt$out, "score_multiz_SARS-CoV-2.txt"),
            quote=F, sep="\t", row.names=F)
