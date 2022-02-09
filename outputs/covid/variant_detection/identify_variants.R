rm(list=ls())

library(here)

variant_detection_dir <- file.path(here(), "outputs", "covid", "variant_detection")
load(file.path(variant_detection_dir, "delta_variants.Rda"))
load(file.path(variant_detection_dir, "omicron_variants.Rda"))

num_delta_variants <- 3052
num_omicron_variants <- 309

# load common variants ----------------------------------------------------

omicron_variants <- unique(omicron_common_variants[, c("type", "pos", "variant", "genome.x")])
omicron_variants$freq <- omicron_variants$genome.x / num_omicron_variants
omicron_variants_substitutions <- subset(omicron_variants, type=="*")
omicron_variants_substitutions$label <- with(omicron_variants_substitutions, paste(pos, type, sep="_"))

delta_variants <- unique(delta_common_variants[, c("type", "pos", "variant", "genome.x")])
delta_variants$freq <- delta_variants$genome.x / num_delta_variants
delta_variants_substitutions <- subset(delta_variants, type=="*")
delta_variants_substitutions$label <- with(delta_variants_substitutions, paste(pos, type, sep="_"))


# identify unique common substitutions ------------------------------------

# n = 29 substitutions unique to delta
delta_unique_substitutions <- subset(delta_variants_substitutions,
                                     !(label %in% omicron_variants_substitutions$label))
delta_unique_substitutions <- subset(delta_unique_substitutions, variant != "n")
delta_unique_substitutions$detect <- "delta"
# n = 47 substitutions unique to omicron
omicron_unique_substitutions <- subset(omicron_variants_substitutions,
                                       !(label %in% delta_variants_substitutions$label))
omicron_unique_substitutions <- subset(omicron_unique_substitutions, variant!="n")
omicron_unique_substitutions$detect <- "omicron"

unique_substitutions <- rbind(delta_unique_substitutions,
                              omicron_unique_substitutions)
unique_substitutions$start <- unique_substitutions$pos - 20

# guide design ------------------------------------------------------------

all_guides <- read.table(file.path(here(), "outputs", "covid", "cas13a_20nt",
                                   "cas13a_results_summary.txt"))
unique_substitutions <- merge(unique_substitutions, all_guides, by="start")

# annotate first position of antitag (n=76)
unique_substitutions$antitag_pos1 <- substr(unique_substitutions$antitag, 1, 1)
# sensitivity >= 95% (n=74)
unique_substitutions <- subset(unique_substitutions, sensitivity_01 >= 0.95)
# specificity = 100% (n=55)
unique_substitutions <- subset(unique_substitutions, specificity == 1)
# no hits to human transcirptome (n=40)
unique_substitutions <- subset(unique_substitutions, match_against_hg38 == 0)
# good crRNA secondary structure (n=14)
unique_substitutions <- subset(unique_substitutions, crRNA_spacer_basepairs == 0)

write.csv(unique_substitutions,
          file=file.path(variant_detection_dir, "unique_substitutions.csv"),
          row.names=F)
