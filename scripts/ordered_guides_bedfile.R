rm(list=ls())

library(here)

# load all guides
plusStrand <- read.table("outputs/cas13a_20nt/cas13a_results_summary.txt")
# minusStrand <- read.table("outputs/cas13a_minusStrand_20nt/cas13a_minusStrand_results_summary.txt")
# dat <- rbind(plusStrand, minusStrand)

# grab ordered guides
inventory <- read.csv("ref_data/NCR inventory - Oligos.csv", stringsAsFactors=F)
inventory <- subset(inventory, Category=="Cas13a crRNA" & grepl("NCR", Oligo.ID) & nchar(Sequence..5..3..)==52)
inventory$spacer <- substr(inventory$Sequence..5..3.., start=(52-20+1), stop=52)
inventory$hairpin <- substr(inventory$Sequence..5..3.., start=1, stop=52-20)

# convert to bed file
bed_file <- data.frame(chr="NC_045512v2",
                       start=plusStrand$start[match(inventory$spacer, plusStrand$spacer)],
                       stop=plusStrand$start[match(inventory$spacer, plusStrand$spacer)]+20-1,
                       strand="+",
                       name=inventory$Oligo.ID)
write.table(bed_file, "outputs/ordered_guides.bed",
            quote=F, row.names=F, col.names=F, sep="\t")
