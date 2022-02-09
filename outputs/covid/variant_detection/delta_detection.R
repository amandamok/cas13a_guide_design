rm(list=ls())

library(here)
library(ggplot2)

min_freq <- 0.75 # frequence to call mutations common to all delta genomes
fill_colors <- RColorBrewer::brewer.pal(3, "Set1")
names(fill_colors) <- c("*", "-", "+")
fill_labels <- c("mismatch", "deletion", "insertion")

parse_mutations <- function(paf) {
  # parse line from minimap paf file
  # output mismatch/indel positions wrt to target
  ## paf: character: line of output from minimap paf file
  paf <- strsplit(paf, split="\t")[[1]]
  align_start <- as.numeric(paf[8]) + 1 # paf is 0-based
  cs_tag <- paf[length(paf)]
  cs_tag <- sub("cs:Z:", "", cs_tag)
  block_seq <- strsplit(cs_tag, split="(=|-|\\+|\\*)")[[1]][-1]
  block_labels <- gsub("(A|a|G|g|T|t|C|c|N|n)", "", cs_tag)
  block_labels <- strsplit(block_labels, split="")[[1]]
  block_lengths <- nchar(block_seq)
  block_lengths[block_labels=="*"] <- 1 # substitution
  block_lengths[block_labels=="+"] <- 0 # insertion (wrt to target)
  mutation_blocks <- which(block_labels != "=")
  if(length(mutation_blocks) == 0) {
    return(NULL)
  } else {
    mutations <- lapply(mutation_blocks,
                        function(x) {
                          tmp_label <- block_labels[x]
                          if(tmp_label %in% c("-", "*")) { # mismatch or deletion
                            # return positions in target corresponding to mismatch/deletion
                            tmp_length <- block_lengths[x] - 1
                            tmp_start <- align_start + sum(block_lengths[1:(x-1)])
                            tmp_positions <- seq(tmp_start, tmp_start + tmp_length)
                            tmp_positions <- data.frame(pos=tmp_positions,
                                                        type=rep(tmp_label, length(tmp_positions)))
                            if(tmp_label == "-") {
                              tmp_positions$variant="-"
                            } else {
                              tmp_positions$variant <- strsplit(block_seq[x], split="")[[1]][2]
                            }
                          } else { # insertion
                            # return positions flanking insertion
                            tmp_start <- align_start + sum(block_lengths[1:(x-1)]) - 1
                            tmp_positions <- tmp_start + 0.5
                            tmp_positions <- data.frame(pos=tmp_positions,
                                                        type=rep(tmp_label, length(tmp_positions)),
                                                        variant=strsplit(block_seq[x], split="")[[1]])
                          }
                          return(tmp_positions)
                        })
    mutations <- do.call(rbind, mutations)
    mutations$genome <- paf[1]
    return(mutations)
  }
}

# load paf output from minimap
delta_paf_fname <- file.path(here(), "ref_data",
                               "gisaid_delta_2022_02_08_09_minimap.paf")
delta_paf <- readLines(delta_paf_fname)

# identify variant names
delta_genomes <- system(paste("cut -f1", delta_paf_fname), intern=T)
delta_genomes <- unique(delta_genomes)
num_genomes <- length(delta_genomes)

# identify mutations relative to reference genome
delta_variants <- lapply(delta_genomes,
                           function(x) {
                             tmp_align <- delta_paf[grepl(x, delta_paf, fixed=T)]
                             tmp_mutations <- lapply(tmp_align, parse_mutations)
                             tmp_mutations <- do.call(rbind, tmp_mutations)
                             return(unique(tmp_mutations))
                           })
delta_variants <- do.call(rbind, delta_variants)

# load bed file of nextstrain genes
genes_bed <- read.table(file.path(here(), "ref_data", "nextstrain_genes.bed"),
                        stringsAsFactors=F,
                        col.names=c("genome", "start", "end", "gene"))
genes_bed$width <- with(genes_bed, end - start)
genes_bed$row <- c(1,2,1,2,1,1,2,1,2,1,2,1)
genes_bed$ymax <- 1.3 - 0.1*(genes_bed$row - 1)
genes_bed$ymin <- 1.21 - 0.1*(genes_bed$row - 1)
genes_bed$width_per_char <- with(genes_bed, width / nchar(gene))

# plot all mutations by type
delta_mutation_summary <- aggregate(genome ~ type + pos,
                                      unique(delta_variants[, c("pos", "type", "genome")]),
                                      length)
summary_plot <- ggplot(delta_mutation_summary,
                       aes(x=pos, y=genome/num_genomes, fill=type)) +
  geom_col(position="stack", width=50) + theme_classic() +
  xlab("position") + ylab("frequency") + expand_limits(y=max(genes_bed$ymax)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  ggtitle(paste0(sum(delta_mutation_summary$genome/num_genomes >= min_freq),
                 " variants present in >", min_freq*100,
                 "% sequenced delta genomes (n =",
                 num_genomes, " )")) +
  scale_fill_manual(values=fill_colors, labels=fill_labels) +
  scale_x_continuous(breaks=seq(0, 30000, 5000)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  geom_hline(yintercept=min_freq, linetype="dotted")
for(x in seq(nrow(genes_bed))) {
  summary_plot <- summary_plot +
    annotate(geom="rect", fill="darkgrey",
             xmin=genes_bed$start[x], xmax=genes_bed$end[x],
             ymin=genes_bed$ymin[x], ymax=genes_bed$ymax[x]) +
    annotate(geom="text", col="black", label=genes_bed$gene[x],
             angle=ifelse(genes_bed$width_per_char[x] <= 200, 90, 0),
             x=with(genes_bed, mean(c(start[x], end[x]))),
             y=with(genes_bed, mean(c(ymin[x], ymax[x]))))
}
summary_plot

# subset to common mutations
common_variants <- subset(delta_mutation_summary, genome/num_genomes >= min_freq)
delta_common_variants <- dplyr::left_join(common_variants, delta_variants,
                                            by=c("type", "pos"))

# load guides
guide_features <- read.csv(file.path(here(), "NCR", "guide_features.csv"), stringsAsFactors=F)
multiplex_32 <- readLines(file.path(here(), "NCR", "32_Guide_List.csv"))
ottlab_8 <- readLines(file.path(here(), "NCR", "Ott8.csv"))[-1]
guides <- subset(guide_features,
                 NCR.id %in% c(paste0("NCR_", multiplex_32), ottlab_8))
guides$width <- nchar(guides$spacer)

# find overlap between NCR guides and common mutations
guide_mutations <- lapply(seq(nrow(guides)),
                          function(x) {
                            tmp_guide <- guides$NCR.id[x]
                            tmp_start <- guides$start[x]
                            tmp_end <- tmp_start + guides$width[x] - 1 + 1 # add first position of anti-tag
                            num_mutated_positions <- sapply(delta_common_variants$pos,
                                                            function(y) {
                                                              sum(y %in% seq(tmp_start, tmp_end))
                                                            })
                            return(num_mutated_positions)
                          })
guide_mutations <- do.call(cbind, guide_mutations)
colnames(guide_mutations) <- guides$NCR.id # no overlap between guides and delta variants

save(delta_variants, delta_common_variants, delta_mutation_summary,
     file=file.path(here(), "outputs", "covid", "variant_detection",
                    "delta_variants.Rda"))
