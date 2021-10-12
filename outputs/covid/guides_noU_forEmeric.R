rm(list=ls())

library(here)

guides <- read.table(file.path(here(), "outputs", "cas13a_20nt", "cas13a_results_summary.txt"), 
                     stringsAsFactors=F)

# choose guides w/o U in spacer
guides <- subset(guides, !grepl("U", spacer))

# remove guides with bad secondary structure
guides <- subset(guides, has_crRNA_hairpin)
guides <- subset(guides, crRNA_spacer_basepairs <= 4)

# remove guides with low specificity
guides <- subset(guides, specificity==1)

# aggregate into overlapping regions
guides$region <- NA
guides$region[1] <- 1
for(x in 2:nrow(guides)) {
  guides$region[x] <- ifelse((guides$start[x] - guides$start[x-1]) <= 20,
                             guides$region[x-1],
                             guides$region[x-1]+1)
}

write.csv(guides, 
          file=file.path(here(), "outputs", "cas13a_20nt", "guides_noU_20201026.csv"))

