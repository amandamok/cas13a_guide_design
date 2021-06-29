library(here)
library(Biostrings)

variant_dir <- file.path(here(), "ref_data", "twist_variants")

wuhcor1 <- readDNAStringSet(file.path(here(), "ref_data", "NC_045512v2.fa"))

guide_features <- read.csv(file.path(here(), "NCR", "guide_features.csv"))
multiplex_32 <- readLines(file.path(here(), "NCR", "32_Guide_List.csv"))
multiplex_32 <- guide_features[match(paste0("NCR_", multiplex_32),
                                     guide_features$NCR.id),]
multiplex_32$target <- sapply(multiplex_32$spacer,
                              function(x) {
                                gsub("U", "T", as.character(reverseComplement(RNAStringSet(x))))
                              })

twist_variant_names <- grep("fasta", list.files(variant_dir), value=T)
twist_variants <- lapply(twist_variant_names,
                         function(x) {
                           readDNAStringSet(file.path(variant_dir, x))
                         })
names(twist_variants) <- sub("\\.fasta", "", twist_variant_names)

twist_variants_alignments <- lapply(twist_variants,
                                    function(x) {
                                      pairwiseAlignment(wuhcor1, x)
                                    })
twist_variants_mismatches <- lapply(twist_variants_alignments, mismatchTable)
twist_variants_indels <- lapply(twist_variants_alignments, indel)

# no mismatches targeted
multiplex_32_mismatches <- t(sapply(multiplex_32$start,
                                    function(x) {
                                      tmp_range <- seq(x, x+20-1)
                                      sapply(twist_variants_mismatches,
                                             function(y) {
                                               sum(y$PatternStart %in% tmp_range)
                                             })
                                    }))
# no indels targeted
multiplex_32_insertion <- t(sapply(multiplex_32$start,
                                   function(x) {
                                     tmp_range <- seq(x, x+20-1)
                                     sapply(twist_variants_indels,
                                            function(y) {
                                              sum(start(insertion(y)) %in% tmp_range)
                                            })
                                   }))
# no deletions targeted
multiplex_32_deletion <- t(sapply(multiplex_32$start,
                                  function(x) {
                                    tmp_range <- seq(x, x+20-1)
                                    sapply(twist_variants_indels,
                                           function(y) {
                                             sum(start(deletion(y)) %in% tmp_range)
                                           })
                                  }))
