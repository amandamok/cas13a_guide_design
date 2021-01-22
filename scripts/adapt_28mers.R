rm(list=ls())

library(here)
library(Biostrings)

ordered <- read.table(file.path(here(), "outputs", "ordered_guides.bed"), 
                      col.names=c("chr", "start", "end", "id"))
ordered <- ordered[order(ordered$start, decreasing=F),]

primary <- read.csv(file.path(here(), "outputs", "NCR_Guide_Data - Data.csv"))

adapt <- subset(ordered, ordered$id %in% primary$NCR.ID)
adapt$start <- adapt$start - 8

genome <- readLines(file.path(here(), "ref_data", "NC_045512v2.fa"))
genome <- genome[!grepl(">", genome)]
genome <- paste(genome, collapse="")

adapt$target <- mapply(substr, start=adapt$start, stop=adapt$end, 
                       MoreArgs=list(x=genome))

write.table(adapt, quote=F, row.names=F, col.names=F, sep="\t",
            file=file.path(here(), "outputs", "adapt_28mers.tsv"))
