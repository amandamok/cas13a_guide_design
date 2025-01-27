rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
ref_dir <- file.path(here(), "ref_data")

# load genomic coordinates for SARS-CoV-2 genes ---------------------------

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

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))

# omit outlier: NCR_1344
guide_rate <- subset(guide_rate, NCR.id != "NCR_1344")

# label top/bottom 10% performing guides
guide_rate_quantiles <- quantile(guide_rate$Estimate, probs=c(0.1, 0.9))
guide_rate$quantile <- "other"
guide_rate$quantile[guide_rate$Estimate <= guide_rate_quantiles["10%"]] <- "< 10%"
guide_rate$quantile[guide_rate$Estimate >= guide_rate_quantiles["90%"]] <- "> 90%"

# generate plot -----------------------------------------------------------

quantile_colors <-setNames(c("blue", "orange", "grey35"),
                           c("< 10%", "> 90%", "other"))

summary_plot <- ggplot(guide_rate, aes(x=start, y=Estimate, col=quantile)) + 
  theme_classic(base_size=8) + geom_point() + 
  geom_errorbar(aes(ymin=Estimate+qnorm(0.025)*Std.Error,
                    ymax=Estimate+qnorm(0.975)*Std.Error)) + 
  xlab("genomic position") + ylab("activator-dependent rate\n(RFU/min)") + 
  coord_cartesian(ylim=c(-11, 105)) + 
  scale_color_manual(values=quantile_colors) + theme(legend.position="none")
gene_plot <- ggplot(gene_bed, aes(xmin=thickStart, xmax=thickEnd, ymin=ymax-4.5, ymax=ymax)) + 
  theme_void() + geom_rect() + 
  geom_text(data=subset(gene_bed, thickEnd-thickStart > 300*nchar(geneName)),
            aes(x=(thickStart + thickEnd)/2, y=ymax-(4.5/2), label=name),
            col="white", size=2.5)
histogram_plot <- ggplot(guide_rate, aes(y=Estimate, fill=quantile)) + 
  theme_classic(base_size=10) + geom_histogram(binwidth=5) + xlab("") + ylab("") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  scale_x_continuous(breaks=c(0, 20)) + coord_cartesian(ylim=c(-25, 90)) + 
  scale_fill_manual(values=quantile_colors)

figure_2A <- gene_plot + summary_plot + plot_spacer() + histogram_plot + 
  plot_layout(byrow=F, widths=c(10,1), heights=c(1,5))

ggsave(filename=file.path(figure_dir, "figure_2A.pdf"),
       plot=figure_2A,
       device="pdf", width=6.5, height=2, units="in")
