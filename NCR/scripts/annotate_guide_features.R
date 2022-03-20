rm(list=ls())

library(here)
library(Biostrings)

project_dir <- file.path(here(), "NCR")
ref_dir <- file.path(here(), "ref_data")

repeat_sequence <- "uagaccaccccaaaaaugaaggggacuaaaac"

# functions ---------------------------------------------------------------

load_sam <- function(sam_fname) {
  # sam_fname: character; file.path to sam alignment file
  sam_colnames <- c("qname", "flag", "rname", "pos", "mapq",
                    "cigar", "rnext", "pnext", "tlen", "seq", "qual")
  sam_data <- read.table(sam_fname, comment.char="@")
  colnames(sam_data)[seq(length(sam_colnames))] <- sam_colnames
  return(sam_data)
}

# guide design pipeline features ------------------------------------------

# 20mers
batch1 <- read.delim(file.path(project_dir, "primary_guide_design_20200710.tsv"))
batch2 <- read.delim(file.path(project_dir, "primary_control_guides_20210129.tsv"))
batch2$spacer <- substr(batch2$sequence, 33, 52)

# 28mers 
adapt_28mer <- read.delim(file.path(project_dir, "adapt_28mers.tsv"))
guide_28mer <- read.delim(file.path(project_dir, "cas13a_28nt_results_summary.txt"))
guide_28mer <- guide_28mer[match(adapt_28mer$target.sequences, guide_28mer$target),]
guide_28mer$NCR.id <- adapt_28mer$New.Name

# merge
guide_features <- intersect(intersect(colnames(batch1), colnames(batch2)),  colnames(guide_28mer))
guide_features <- rbind(batch1[, guide_features], batch2[, guide_features], guide_28mer[, guide_features])

# crRNA secondary structure -----------------------------------------------

guide_features$spacer_structure <- guide_features$crRNA_spacer_basepairs / nchar(guide_features$spacer)
guide_features$crRNA_spacer_basepairs <- factor(guide_features$crRNA_spacer_basepairs, 
                                                levels=sort(unique(guide_features$crRNA_spacer_basepairs)))

# bin guides by crRNA secondary structure 
# has_hairpin & crRNA_spacer_basepairs
guide_features$group <- sapply(seq(nrow(guide_features)),
                               function(x) {
                                 if(guide_features$has_crRNA_hairpin[x] & 
                                    guide_features$crRNA_spacer_basepairs[x] == 0) {
                                   return("good guide")
                                 } else {
                                   paste0(ifelse(guide_features$has_crRNA_hairpin[x], 
                                                 "hasHairpin", "noHairpin"), "/",
                                          guide_features$crRNA_spacer_basepairs[x], "bp")
                                 }
                               })
guide_features$group <- factor(guide_features$group,
                               levels=c("good guide", 
                                        "hasHairpin/4bp", "noHairpin/4bp",
                                        "hasHairpin/6bp", "noHairpin/7bp",
                                        "noHairpin/8bp", "noHairpin/9bp",
                                        "hasHairpin/10bp", "noHairpin/10bp",
                                        "noHairpin/11bp", "hasHairpin/12bp",
                                        "noHairpin/13bp", "noHairpin/14bp",
                                        "noHairpin/15bp", "hasHairpin/16bp", 
                                        "noHairpin/16bp","noHairpin/17bp"))
guide_features$class <- sapply(as.numeric(sub("NCR_", "", guide_features$NCR.id)),
                               function(x) {
                                 if(x >= 504 & x <= 599) {
                                   return("good")
                                 } else {
                                   if(x >= 600 & x <= 614) {
                                     return("8plex")
                                   } else {
                                     if(x >= 1305 & x <= 1352) {
                                       return("bad")
                                     } else {
                                       if(x >= 1353 & x <= 1400) {
                                         return("random")
                                       } else {
                                         if(x >= 1407 & x <= 1418) {
                                           return("28mer")
                                         } else {
                                           return(NA)
                                         }
                                       }
                                     }
                                   }
                                 }
                               })
guide_features$class <- factor(guide_features$class, levels=c("good", "8plex", "bad", "random", "28mer"))

# load crRNA MFE from RNAfold
all_guides <- read.table(file.path(here(), "outputs", "covid", "cas13a_20nt",
                                   "cas13a_results_summary.txt"), header=T)
guide_features$gRNA_MFE <- all_guides$gRNA_MFE[match(guide_features$start,
                                                     all_guides$start)]

# load crRNA dot bracket structure
RNAfold_20nt <- read.table(file.path(here(), "outputs", "covid", "cas13a_20nt", 
                                     "score_RNAfold_crRNAs.txt"), header=T)
RNAfold_28nt <- read.table(file.path(here(), "outputs", "covid", "cas13a_28nt",
                                     "score_RNAfold_crRNAs.txt"), header=T)
guide_features$structure <- sapply(seq(nrow(guide_features)),
                                   function(x) {
                                     tmp_start <- guide_features$start[x]
                                     tmp_spacer_length <- nchar(guide_features$spacer[x])
                                     tmp_structure <- subset(get(paste0("RNAfold_", 
                                                                        tmp_spacer_length, 
                                                                        "nt")),
                                                             start == tmp_start)
                                     return(tmp_structure$structure[1])
                                   })
guide_features$direct_repeat <- factor(substr(guide_features$structure, 
                                              1, nchar(repeat_sequence)))
guide_features$direct_repeat <- relevel(guide_features$direct_repeat,
                                        ref=".......((((.........))))........")

# load target-crRNA hybridization MFE
guide_features$target <- as.character(reverseComplement(RNAStringSet(guide_features$spacer)))
guide_features$hybridization_MFE <- sapply(seq(nrow(guide_features)),
                                           function(x) {
                                             tmp <- system(paste0('printf "', 
                                                                  guide_features$spacer[x],
                                                                  '\n',
                                                                  guide_features$target[x],
                                                                  '" | RNAduplex'),
                                                           intern=T)
                                             tmp <- strsplit(tmp, split=" ")[[1]][[11]]
                                           })
guide_features$hybridization_MFE <- as.numeric(sub("\\(", "",
                                                   sub("\\)", "", 
                                                       guide_features$hybridization_MFE)))

# viral target secondary structure ----------------------------------------

viral_genome <- readDNAStringSet(file.path(ref_dir, "NC_045512v2.fa"))
viral_genome <- data.frame(pos = seq(nchar(viral_genome)),
                           base = unlist(strsplit(as.character(viral_genome), 
                                                  split="")),
                           spacer_base = unlist(strsplit(as.character(complement(viral_genome)), 
                                                         split="")),
                           row.names=NULL)
viral_genome$base[viral_genome$base == "T"] <- "U"
viral_genome$spacer_base[viral_genome$spacer_base == "T"] <- "U"

# Manfredonia et al. NAR (2020)
manfredonia_prefix <- "manfredonia_2020_lowShannon_highSHAPE"
manfredonia_sam_fname <- file.path(ref_dir, paste0(manfredonia_prefix, ".sam"))
if(!file.exists(manfredonia_sam_fname)) {
  manfredonia <- openxlsx::read.xlsx(file.path(ref_dir, paste0(manfredonia_prefix, ".xlsx")),
                                     colNames=T, startRow=2)
  manfredonia_fasta <- Biostrings::DNAStringSet(manfredonia$Sequence)
  names(manfredonia_fasta) <- with(manfredonia, paste("region", From, To, sep="_"))
  manfredonia_fasta_fname <- file.path(ref_dir, paste0(manfredonia_prefix, ".fa"))
  Biostrings::writeXStringSet(manfredonia_fasta, filepath=manfredonia_fasta_fname)
  system(paste("bowtie --norc -v 2 -S -f", file.path(ref_dir, "wuhCor1"),
               manfredonia_fasta_fname, ">", manfredonia_sam_fname))
}
## add annotation
manfredonia_sam <- load_sam(manfredonia_sam_fname)
viral_genome$manfredonia <- "."
for(x in seq(nrow(manfredonia_sam))) {
  region_start <- manfredonia_sam$pos[x]
  region_end <- region_start + nchar(manfredonia_sam$seq[x]) - 1
  viral_genome$manfredonia[region_start:region_end] <- "lowShannon/highSHAPE"
}

# Lan et al. bioRxiv (2020)
lan_structured_prefix <- "lan_2020_structured"
lan_structured_sam_fname <- file.path(ref_dir, paste0(lan_structured_prefix, ".sam"))
if(!file.exists(lan_structured_sam_fname)) {
  lan_structured <- read.csv(file.path(ref_dir, paste0(lan_structured_prefix, ".csv")))
  lan_structured <- cbind("MN985325.1", lan_structured)
  lan_structured_bed_fname <- file.path(ref_dir, paste0(lan_structured_prefix, ".bed"))
  write.table(lan_structured, file=lan_structured_bed_fname,
              sep="\t", quote=F, row.names=F, col.names=F)
  lan_structured_fasta_fname <- file.path(ref_dir, paste0(lan_structured_prefix, ".fa"))
  system(paste("bedtools getfasta -fi",
               file.path(ref_dir, "MN985325v1.fa"),
               "-bed", lan_structured_bed_fname,
               "-fo", lan_structured_fasta_fname))
  system(paste("bowtie --norc -v 2 -S -f", file.path(ref_dir, "wuhCor1"),
               lan_structured_fasta_fname, ">", lan_structured_sam_fname))
}
lan_unstructured_prefix <- "lan_2020_unstructured"
lan_unstructured_sam_fname <- file.path(ref_dir, paste0(lan_unstructured_prefix, ".sam"))
if(!file.exists(lan_unstructured_sam_fname)) {
  lan_unstructured <- read.csv(file.path(ref_dir, paste0(lan_unstructured_prefix, ".csv")))
  lan_unstructured <- cbind("MN985325.1", lan_unstructured)
  lan_unstructured_bed_fname <- file.path(ref_dir, paste0(lan_unstructured_prefix, ".bed"))
  write.table(lan_unstructured, file=lan_unstructured_bed_fname,
              sep="\t", quote=F, row.names=F, col.names=F)
  lan_unstructured_fasta_fname <- file.path(ref_dir, paste0(lan_unstructured_prefix, ".fa"))
  system(paste("bedtools getfasta -fi",
               file.path(ref_dir, "MN985325v1.fa"),
               "-bed", lan_unstructured_bed_fname,
               "-fo", lan_unstructured_fasta_fname))
  system(paste("bowtie --norc -v 2 -S -f", file.path(ref_dir, "wuhCor1"),
               lan_unstructured_fasta_fname, ">", lan_unstructured_sam_fname))
}
# add annotations
lan_structured_sam <- load_sam(lan_structured_sam_fname)
lan_unstructured_sam <- load_sam(lan_unstructured_sam_fname)
viral_genome$lan <- "."
for(x in seq(nrow(lan_structured_sam))) {
  region_start <- lan_structured_sam$pos[x]
  region_end <- region_start + nchar(lan_structured_sam$seq[x]) - 1
  viral_genome$lan[region_start:region_end] <- "structured"
}
for(x in seq(nrow(lan_unstructured_sam))) {
  region_start <- lan_unstructured_sam$pos[x]
  region_end <- region_start + nchar(lan_unstructured_sam$seq[x]) - 1
  region <- viral_genome$lan[region_start:region_end]
  region <- paste0(region, "_unstructured")
  viral_genome$lan[region_start:region_end] <- region
}
viral_genome$lan[viral_genome$lan=="._unstructured"] <- "unstructured"
viral_genome$lan[viral_genome$lan=="structured_unstructured"] <- "both"

# Sun et al. Cell (2021)
sun_invivo <- openxlsx::read.xlsx(file.path(ref_dir, "sun_2021_icSHAPE.xlsx"),
                                  sheet="SARS2-invivo", na.strings="NULL")
viral_genome$sun_invivo <- sun_invivo$`icSHAPE-score`
sun_invitro <- openxlsx::read.xlsx(file.path(ref_dir, "sun_2021_icSHAPE.xlsx"),
                                   sheet="SARS2-invitro", na.strings="NULL")
viral_genome$sun_invitro <- sun_invitro$`icSHAPE-score`

# Huston et al. Molecular Cell*(2021)
huston <- read.table(file.path(ref_dir, "huston_2021_structure.ct"), skip=1,
                     col.names=c("index", "base", "previous_index", "next_index", 
                                 "paired_to", "natural_numbering"))
viral_genome$huston <- as.numeric(huston$paired_to)

# clustering of in vivo metrics
viral_genome_invivo <- data.frame(lan = (viral_genome$lan %in% c("both", "unstructured")),
                                  sun = (viral_genome$sun_invivo >= 0.5),
                                  huston = (viral_genome$huston == 0))
viral_genome$invivo <- ifelse(rowSums(viral_genome_invivo, na.rm=T) >= 2,
                              "unstructured", "structured")

# generate crRNA annotations
guide_structure <- sapply(seq(nrow(guide_features)),
                          function(x) {
                            region_start <- guide_features$start[x]
                            region_end <- region_start + nchar(guide_features$spacer[x]) - 1
                            region_features <- viral_genome[match(region_start:region_end,
                                                                  viral_genome$pos),]
                            region_annotations <- c(manfredonia = mean(region_features$manfredonia == "lowShannon/highSHAPE"),
                                                    lan = mean(region_features$lan == "unstructured"),
                                                    sun_invivo = mean(region_features$sun_invivo),
                                                    sun_invitro = mean(region_features$sun_invitro),
                                                    huston = mean(region_features$huston == 0),
                                                    invivo = mean(region_features$invivo == "unstructured"))
                            return(region_annotations)
                          })
guide_structure <- data.frame(t(guide_structure))
guide_structure$invivo_cat <- sapply(guide_structure$invivo,
                                     function(x) {
                                       ifelse(x==0,
                                              "DS",
                                              ifelse(x==1,
                                                     "SS",
                                                     "SS-DS"))
                                     })
guide_features <- cbind(guide_features, guide_structure)

# neighborhood of target region -------------------------------------------

U_window <- 30
guide_features$downstream_U <- sapply(seq(nrow(guide_features)),
                                      function(x) { # downstream wrt spacer
                                        # region_start <- guide_features$start[x] +
                                        #   nchar(guide_features$spacer[x])
                                        # region_end <- region_start + U_window - 1
                                        region_end <- guide_features$start[x] - 1
                                        region_start <- region_end - U_window + 1
                                        tmp_region <- viral_genome[region_start:region_end,]
                                        mean(tmp_region$base == "U")
                                      })
guide_features$downstream_unstructured_U <- sapply(seq(nrow(guide_features)),
                                                   function(x) { # downstream wrt spacer
                                                     # region_start <- guide_features$start[x] +
                                                     #   nchar(guide_features$spacer[x])
                                                     # region_end <- region_start + U_window - 1
                                                     region_end <- guide_features$start[x] - 1
                                                     region_start <- region_end - U_window + 1
                                                     tmp_region <- viral_genome[region_start:region_end,]
                                                     mean((tmp_region$base == "U") &
                                                            (tmp_region$invivo == "unstructured"))
                                                   })
guide_features$U_context <- sapply(seq(nrow(guide_features)),
                                   function(x) { # upstream/downstream wrt spacer
                                     upstream_start <- guide_features$start[x] - 1
                                     upstream_end <- upstream_start - U_window + 1
                                     downstream_start <- guide_features$start[x] + 
                                       nchar(guide_features$spacer[x])
                                     downstream_end <- downstream_start + U_window - 1
                                     tmp_region <- rbind(viral_genome[upstream_start:upstream_end,],
                                                         viral_genome[downstream_start:downstream_end,])
                                     mean(tmp_region$base=="U")
                                   })
guide_features$U_context_unstructured <- sapply(seq(nrow(guide_features)),
                                                function(x) {
                                                  upstream_start <- guide_features$start[x] - 1
                                                  upstream_end <- upstream_start - U_window + 1
                                                  downstream_start <- guide_features$start[x] + 
                                                    nchar(guide_features$spacer[x])
                                                  downstream_end <- downstream_start + U_window - 1
                                                  tmp_region <- rbind(viral_genome[upstream_start:upstream_end,],
                                                                      viral_genome[downstream_start:downstream_end,])
                                                  mean((tmp_region$base=="U") &
                                                         (tmp_region$invivo=="unstructured"))
                                                })
guide_features$antitag <- sapply(seq(nrow(guide_features)),
                                 function(x) {
                                   tmp_position <- guide_features$start[x] + 
                                     nchar(guide_features$spacer[x])
                                   paste(viral_genome$base[tmp_position+(0:3)],
                                         collapse="")
                                 })
guide_features$antitag_pos1 <- substr(guide_features$antitag, 1, 1)
guide_features$antitag_label <- sapply(guide_features$antitag,
                                       function(x) {
                                         ifelse(grepl("^GUU", x), "GUU",
                                                ifelse(grepl("^GU", x), "GU",
                                                       ifelse(grepl("^G", x), "G",
                                                              substr(x, 1, 1))))
                                       })
guide_features$antitag_label <- factor(guide_features$antitag_label,
                                       levels=c("A", "C", "U", "G", "GU", "GUU"))

# cis cleavage site -------------------------------------------------------

cis_cleavage_start <- 20
cis_cleavage_window <- 20
guide_features$downstream_U <- sapply(seq(nrow(guide_features)),
                                      function(x) {
                                        region_start <- guide_features$start[x] +
                                          nchar(guide_features$spacer[x])
                                        region_end <- region_start + U_window - 1
                                        tmp_region <- viral_genome[region_start:region_end,]
                                        mean(tmp_region$base == "U")
                                      })
guide_features$cis_cleavage_site <- sapply(seq(nrow(guide_features)),
                                           function(x) {
                                             region_end <- guide_features$start[x] - 1 - 
                                               cis_cleavage_start + cis_cleavage_window/2
                                             region_start<- region_end - cis_cleavage_window + 1
                                             tmp_region <- viral_genome[region_start:region_end,]
                                             paste(tmp_region$base, collapse="")
                                           })

# write annotations to file -----------------------------------------------

write.csv(guide_features, 
          file=file.path(project_dir, "data", "guide_features.csv"), 
          row.names=F)
write.csv(viral_genome, 
          file=file.path(project_dir, "data", "viral_genome_features.csv"), 
          row.names=F)
