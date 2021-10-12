rm(list=ls())

library(here)
library(Biostrings)
library(ggplot2)
library(patchwork)

twist_variants_dir <- file.path(here(), "ref_data", "twist_variants")
delta_variants <- c("EPI_ISL_1544014", "EPI_ISL_2695467", "EPI_ISL_2693246")

# load sequences
delta_variants_seq <- lapply(delta_variants,
                             function(x) {
                               readDNAStringSet(file.path(twist_variants_dir,
                                                          paste0(x, ".fasta")))
                             })
delta_variants_seq <- DNAStringSet(sapply(delta_variants_seq, as.character))
wildtype <- readDNAStringSet(file.path(twist_variants_dir, 
                                       "MN908947.3.fasta"))

# find mutations between delta and wildtype
delta_alignments <- pairwiseAlignment(delta_variants_seq, wildtype)
delta_alignments <- c(as.character(wildtype),
                      as.character(alignedPattern(delta_alignments)))
delta_alignments <- data.frame(seq(nchar(wildtype)),
                               matrix(unlist(strsplit(delta_alignments,
                                                      split="")),
                                      ncol=4, byrow=F))
colnames(delta_alignments) <- c("position", 
                                strsplit(names(wildtype), split=" ")[[1]][1],
                                delta_variants)
for(x in delta_variants) {
  delta_alignments[paste0(x, "_variant")] <- (delta_alignments$MN908947.3 != delta_alignments[x])
}
delta_alignments$mutation <- rowSums(delta_alignments[, grepl("_variant", colnames(delta_alignments))])
delta_alignments <- subset(delta_alignments, mutation != 0)
delta_alignments <- subset(delta_alignments, 
                           (position >= 100) & 
                             (position <= (nchar(wildtype) - 100)))

# find guides that overlap mutations
all_guides <- read.csv(file.path(here(), "NCR", "guide_features.csv"))
delta_alignments$guides <- sapply(delta_alignments$position,
                                  function(x) {
                                    tmp_guides <- which((all_guides$start <= x) &
                                                          ((all_guides$start + 19) >= x))
                                    tmp_guides <- all_guides$NCR.id[tmp_guides]
                                    if(length(tmp_guides) == 0) {
                                      return("")
                                    } else {
                                      return(paste(tmp_guides, collapse="|"))
                                    }
                                  })
delta_guides <- unique(delta_alignments$guides)[-1]
delta_guides <- unlist(strsplit(delta_guides, split="\\|"))

# plot
guide_rate <- read.csv(file.path(here(), "NCR", "guide_rate.csv"))
guide_rate$num_delta_variants <- delta_alignments$mutation[match(guide_rate$NCR.id, delta_alignments$guides)]
guide_rate$num_delta_variants[guide_rate$NCR.id %in% c("NCR_1312", "NCR_1399")] <- 1
guide_rate$num_delta_variants[is.na(guide_rate$num_delta_variants)] <- 0
guide_rate$delta_alpha <- ifelse(is.na(guide_rate$num_delta_variants), 0.1, 1)
gene_bed <- read.table(file.path(here(), "ref_data", "cov2_genes.txt"),
                       comment.char="@", header=T)
gene_bed$overlap <- c(F, 
                      sapply(2:nrow(gene_bed),
                             function(x) {
                               gene_bed$thickStart[x] < (gene_bed$thickEnd[x-1]+50)
                             }))
gene_bed$ymax <- NA
for(x in seq(nrow(gene_bed))) {
  if(x==1) {
    gene_bed$ymax[x] <- -15
  } else {
    gene_bed$ymax[x] <- ifelse(gene_bed$overlap[x],
                               ifelse(gene_bed$ymax[x-1]==-15, -20, -15),
                               ifelse(gene_bed$ymax[x-1]==-15, -15, -20))
  }
}

(ggplot(guide_rate, aes(x=start, y=Estimate, col=as.factor(num_delta_variants))) + 
    theme_classic() + geom_point(aes(alpha=delta_alpha)) + 
    geom_errorbar(aes(ymin=Estimate+qnorm(0.05)*Std.Error,
                      ymax=Estimate-qnorm(0.05)*Std.Error), alpha=0.25) + 
    geom_text(data=subset(guide_rate, num_delta_variants > 0), 
              aes(label=NCR.id, y=Estimate + 4), size=3, show_guide=F) + 
    xlab("position along genome") + ylab("activator-dependent rate\n(RFU/min)") + 
    ylim(-10, 100) + guides(alpha="none") + labs(col="# delta variants targeted:") + 
    scale_color_manual(values=c("black", "blue", "red")) + 
    theme(axis.title=element_text(size=10), axis.text=element_text(size=10),
          legend.position="top")) +
  (ggplot(gene_bed, aes(xmin=thickStart, xmax=thickEnd, ymin=ymax-4.5, ymax=ymax)) + 
     theme_void() + geom_rect() + 
     geom_text(data=subset(gene_bed, thickEnd-thickStart > 300*nchar(geneName)),
               aes(x=(thickStart + thickEnd)/2, y=ymax-(4.5/2), label=name),
               col="white", size=2.5)) +
  (ggplot(guide_rate, aes(y=Estimate)) + theme_classic() + 
     geom_histogram(binwidth=5) + xlab("") + ylab("") + ylim(-10, 100) + 
     theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
           axis.text.x=element_text(size=10)) + 
     scale_x_continuous(breaks=c(0, 20))) + 
  plot_spacer() + 
  plot_layout(byrow=F, widths=c(10,1), heights=c(5,1))

write.csv(delta_alignments, file.path(here(), "outputs", "delta_guides.csv"),
          quote=F, row.names=F)
