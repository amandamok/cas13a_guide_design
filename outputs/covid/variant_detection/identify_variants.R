rm(list=ls())

library(here)

# load variant fasta files ------------------------------------------------

variant_dir <- file.path(here(), "ref_data", "twist_variants")
twist_variants <- read.table(file.path(variant_dir, "twist_variants.tsv"), 
                             sep="\t", header=T, stringsAsFactors=F)

wt_control <- subset(twist_variants, Variant=="Control")$GenBank.ID
delta_variants <- subset(twist_variants, WHO.label=="delta")$GenBank.ID
omicron_variants <- subset(twist_variants, WHO.label=="omicron")$GenBank.ID

# generate minimap alignment ----------------------------------------------

for(x in c(delta_variants, omicron_variants)) {
  if(!file.exists(file.path(variant_dir, paste0(x, ".vars")))) {
    paste("minimap2 --cs", 
          file.path(variant_dir, "MN908947.3.fasta"), 
          file.path(variant_dir, paste0(x, ".fasta")), 
          "| paftools.js call -L 10000 - >",
          file.path(variant_dir, paste0(x, ".vars")))
  }
}

# load nextstrain clade definitions ---------------------------------------

# https://github.com/nextstrain/ncov/blob/master/defaults/clades.tsv
clade_snps <- read.table(file.path(here(), "outputs", "covid",
                                   "variant_detection", "clades.tsv"),
                         sep="\t", header=T)
clade_snps$variant <- with(clade_snps, paste(site, alt, sep="_"))

# identify clades ---------------------------------------------------------

detection_clades <- c("Delta", "Omicron")
detection_clades <- grep(paste0("(", paste(detection_clades, collapse=")|("), ")"),
                         unique(clade_snps$clade), value=T)

# identify differentiating snps -------------------------------------------

detection_snps <- subset(clade_snps, clade %in% detection_clades)
detection_snps <- aggregate(clade~variant, detection_snps, unique)
detection_snps <- cbind(detection_snps,
                        sapply(detection_clades,
                               function(x) {
                                 sapply(detection_snps$clade,
                                        function(y) {
                                          as.numeric(x %in% y)
                                        })
                               }))
detection_snps <- subset(detection_snps, select=-clade)
detection_snps$position <- as.numeric(sub("_.", "", detection_snps$variant))
detection_snps$alt_allele <- sub(".*_", "", detection_snps$variant)

detection_snps$num_delta_nextstrain <- rowSums(detection_snps[, grepl("Delta", colnames(detection_snps))])
detection_snps$num_omicron_nextstrain <- rowSums(detection_snps[, grepl("Omicron", colnames(detection_snps))])

# detection_snps$clade_specific <- rowSums(subset(detection_snps,
#                                                 select=grepl("2", colnames(detection_snps)))) == 1
# detection_snps$delta_specific <- (rowSums(subset(detection_snps,
#                                                  select=grepl("Delta", colnames(detection_snps)))) == 
#                                     sum(grepl("Delta", detection_clades))) & 
#   (rowSums(subset(detection_snps,
#                   select=grepl("Omicron", colnames(detection_snps)))) == 0)
# detection_snps$omicron_specific <- (rowSums(subset(detection_snps,
#                                                    select=grepl("Omicron", colnames(detection_snps)))) == 
#                                       sum(grepl("Omicron", detection_clades))) & 
#   (rowSums(subset(detection_snps,
#                   select=grepl("Delta", colnames(detection_snps)))) == 0)
# detection_snps$detect_variant <- with(detection_snps, delta_specific | omicron_specific)

# identify crRNA spacers --------------------------------------------------

all_crRNAs <- read.table(file.path(here(), "outputs", "covid", "cas13a_20nt", 
                                   "cas13a_results_summary.txt"))
all_crRNAs$antitag_pos1 <- substr(all_crRNAs$antitag, 1, 1)

detection_crRNAs <- all_crRNAs[match(detection_snps$position-20,
                                     all_crRNAs$start),]
detection_crRNAs <- cbind(detection_crRNAs, detection_snps)

# check delta variants ----------------------------------------------------

delta_mutations <- data.frame(delta_allele=sapply(seq(nrow(detection_crRNAs)),
                                                  function(x) {
                                                    if(detection_crRNAs$num_delta_nextstrain[x] == 3) {
                                                      # SNP in all 3 Delta clades in Nextstrain
                                                      # return alternative allele
                                                      return(detection_crRNAs$alt_allele[x])
                                                    } else {
                                                      # return WT allele
                                                      return(detection_crRNAs$antitag_pos1[x])
                                                    }
                                                  }))

for(x in delta_variants) {
  tmp_vars <- read.table(file.path(variant_dir, paste0(x, ".vars")),
                         sep="\t", header=F, skip=1)
  colnames(tmp_vars) <- c("type", "chr", "start", "end", "query_depth", 
                          "mapping_quality", "ref_allele", "alt_allele",
                          "query_name", "query_start", "query_end", 
                          "query_orientation")
  tmp_vars$ref_allele[tmp_vars$ref_allele=="t"] <- "u"
  tmp_vars$alt_allele[tmp_vars$alt_allele=="t"] <- "u"
  delta_mutations[x] <- sapply(seq(nrow(detection_crRNAs)),
                               function(y) {
                                 which_row <- match(detection_crRNAs$position[y], 
                                                    tmp_vars$end)
                                 if(is.na(which_row)) {
                                   # SNP not detected
                                   return(detection_crRNAs$antitag_pos1[y])
                                 } else {
                                   return(toupper(tmp_vars$alt_allele[which_row]))
                                 }
                               })
}

delta_mutations$num_delta_twist <- sapply(seq(nrow(delta_mutations)),
                                          function(x) {
                                            sum(delta_mutations$delta_allele[x] == 
                                                  delta_mutations[x, delta_variants])
                                          })

detection_crRNAs <- cbind(detection_crRNAs, delta_mutations)

# check omicron variants --------------------------------------------------

omicron_mutations <- data.frame(omicron_allele=sapply(seq(nrow(detection_crRNAs)),
                                                      function(x) {
                                                        if(detection_crRNAs$num_omicron_nextstrain[x] == 3) {
                                                          # SNP in all 3 Omicron clades in Nextstrain
                                                          # return alternative allele
                                                          return(detection_crRNAs$alt_allele[x])
                                                        } else {
                                                          # return WT allele
                                                          return(detection_crRNAs$antitag_pos1[x])
                                                        }
                                                      }))

for(x in omicron_variants) {
  tmp_vars <- read.table(file.path(variant_dir, paste0(x, ".vars")),
                         sep="\t", header=F, skip=1)
  colnames(tmp_vars) <- c("type", "chr", "start", "end", "query_depth", 
                          "mapping_quality", "ref_allele", "alt_allele",
                          "query_name", "query_start", "query_end", 
                          "query_orientation")
  tmp_vars$ref_allele[tmp_vars$ref_allele=="t"] <- "u"
  tmp_vars$alt_allele[tmp_vars$alt_allele=="t"] <- "u"
  omicron_mutations[x] <- sapply(seq(nrow(detection_crRNAs)),
                                 function(y) {
                                   which_row <- match(detection_crRNAs$position[y], 
                                                      tmp_vars$end)
                                   if(is.na(which_row)) {
                                     # SNP not detected
                                     return(detection_crRNAs$antitag_pos1[y])
                                   } else {
                                     return(toupper(tmp_vars$alt_allele[which_row]))
                                   }
                                 })
}

omicron_mutations$num_omicron_twist <- sapply(seq(nrow(omicron_mutations)),
                                              function(x) {
                                                sum(omicron_mutations$omicron_allele[x] == 
                                                      omicron_mutations[x, omicron_variants])
                                              })

detection_crRNAs <- cbind(detection_crRNAs, omicron_mutations)

# write results -----------------------------------------------------------

detection_crRNAs$good_SNP <- with(detection_crRNAs,
                                  (num_delta_twist==3) & (num_omicron_twist==3))

write.csv(detection_crRNAs,
          file=file.path(here(), "outputs", "covid",
                         "variant_detection", "detection_crRNAs.csv"),
          row.names=F)
