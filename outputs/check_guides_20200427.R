library(here)

bowtie_path <- "/mnt/ingolialab/linux-x86_64/bin/bowtie"

setwd(file.path(here(), "outputs"))

guides <- readLines("guides_20200427.txt")
spacers <- sapply(guides,
                  function(x) {
                    substr(x, start=nchar(x)-20+1, stop=nchar(x))
                  })
spacers <- gsub("U", "T", spacers)
writeLines(spacers, con="spacers_20200427.txt")

spacers_align <- system(paste(bowtie_path, "-S -r ~/covid-19/ref_data/wuhCor1 spacers_20200427.txt"), intern=T)
spacers_align <- spacers_align[!grepl("@", spacers_align)]
spacers_align <- data.frame(matrix(unlist(strsplit(spacers_align, split="\t")), 
                                   nrow=length(spacers_align), byrow=T), 
                            stringsAsFactors=F)

# summary(spacers_align[,10] == spacers)
summary(spacers == as.character(reverseComplement(DNAStringSet(spacers_align[,10]))))