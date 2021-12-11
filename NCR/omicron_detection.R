rm(list=ls())

library(here)
library(ggplot2)

min_freq <- 0.75 # frequence to call mutations common to all omicron genomes
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
  mutation_positions <- lapply(mutation_blocks,
                               function(x) {
                                 tmp_label <- block_labels[x]
                                 if(tmp_label %in% c("-", "*")) { # mismatch or deletion
                                   # return positions in target corresponding to mismatch/deletion
                                   tmp_length <- block_lengths[x] - 1
                                   tmp_start <- align_start + sum(block_lengths[1:(x-1)])
                                   tmp_positions <- seq(tmp_start, tmp_start + tmp_length)
                                 } else { # insertion
                                   # return positions flanking insertion
                                   tmp_start <- align_start + sum(block_lengths[1:(x-1)]) - 1
                                   tmp_end <- tmp_start + 1
                                   tmp_positions <- c(tmp_start, tmp_end)
                                 }
                                 names(tmp_positions) <- rep(tmp_label, length(tmp_positions))
                                 return(tmp_positions)
                               })
  mutation_positions <- unlist(mutation_positions)
}

# load bed file of nextstrain genes
genes_bed <- read.table(file.path(here(), "ref_data", "nextstrain_genes.bed"),
                        stringsAsFactors=F,
                        col.names=c("genome", "start", "end", "gene"))
genes_bed$width <- with(genes_bed, end - start)
genes_bed$row <- c(1,2,1,2,1,1,2,1,2,1,2,1)
genes_bed$ymax <- 1.3 - 0.1*(genes_bed$row - 1)
genes_bed$ymin <- 1.21 - 0.1*(genes_bed$row - 1)
genes_bed$width_per_char <- with(genes_bed, width / nchar(gene))

# load paf output from minimap
omicron_paf_fname <- file.path(here(), "ref_data",
                               "gisaid_omicron_variants_minimap.paf")
omicron_paf <- readLines(omicron_paf_fname)

# identify variant names
omicron_variants <- system(paste("cut -f1", omicron_paf_fname), intern=T)
omicron_variants <- unique(omicron_variants)
num_variants <- length(omicron_variants)

# identify mutations relative to reference genome
omicron_mutations <- lapply(omicron_variants,
                            function(x) {
                              tmp_align <- omicron_paf[grepl(x, omicron_paf, fixed=T)]
                              tmp_mutations <- lapply(tmp_align, parse_mutations)
                              tmp_mutations <- unlist(tmp_mutations)
                              tmp_mutations <- data.frame(variant=x,
                                                          type=names(tmp_mutations),
                                                          position=tmp_mutations)
                              return(unique(tmp_mutations))
                            })
names(omicron_mutations) <- omicron_variants

# plot all mutations by type
mutation_freq <- aggregate(variant ~ position,
                           do.call(rbind, omicron_mutations),
                           length)
mutation_summary <- aggregate(variant ~ type + position,
                              do.call(rbind, omicron_mutations),
                              length)
summary_plot <- ggplot(mutation_summary,
                       aes(x=position, y=variant/num_variants, fill=type)) +
  geom_col(position="stack", width=50) + theme_classic() +
  xlab("position") + ylab("frequency") + expand_limits(y=max(genes_bed$ymax)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  ggtitle(paste0(sum(mutation_freq$variant/num_variants >= min_freq),
                 " variants present in >", min_freq*100,
                 "% sequenced omicron genomes (n =",
                 length(omicron_variants), ")")) +
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
common_mutations <- subset(mutation_freq, variant/num_variants >= min_freq)
omicron_mutations_common <- lapply(omicron_mutations,
                                   function(x) {
                                     subset(x, position %in% common_mutations$position)
                                   })

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
                            num_mutated_positions <- sapply(omicron_mutations_common,
                                                            function(y) {
                                                              sum(y$position %in% seq(tmp_start, tmp_end))
                                                            })
                            return(num_mutated_positions)
                          })
guide_mutations <- do.call(cbind, guide_mutations)
colnames(guide_mutations) <- guides$NCR.id

# compute number of missed targets per guide
non_targeting_guides <- colSums(guide_mutations > 0)
non_targeting_guides <- non_targeting_guides[non_targeting_guides > 0]
non_targeting_guides # only NCR_576

NCR_576_start <- guides$start[guides$NCR.id=="NCR_576"]
NCR_576_end <- guides$start[guides$NCR.id=="NCR_576"] + guides$width[guides$NCR.id=="NCR_576"] - 1
ggplot(subset(mutation_summary,
              position %in% seq(NCR_576_start, NCR_576_end)),
       aes(x=position, y=variant, fill=type)) +
  geom_col(position="stack") + theme_classic() +
  scale_fill_manual(values=fill_colors, labels=fill_labels) +
  xlab("position") + ylab("# mutations") +
  xlim(NCR_576_start-0.5, NCR_576_end+0.5) + ylim(-15, num_variants) +
  annotate(geom="text", x=seq(NCR_576_start, NCR_576_end), y=-12,
           label=strsplit(guides$target[guides$NCR.id=="NCR_576"], split="")[[1]]) +
  theme(legend.position="top")
