rm(list=ls())

library(here)
library(ggplot2)
library(ggpubr)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))
guide_rate$antitag_pos1 <- factor(guide_rate$antitag_pos1,
                                  levels=c("G", "A", "C", "U"))

# compute pvalues ---------------------------------------------------------

antitag_pos1_pvalue <- compare_means(Estimate ~ antitag_pos1,
                                     guide_rate, method="wilcox.test", 
                                     p.adjust.method="none")

# generate plot -----------------------------------------------------------

fill_colors <- RColorBrewer::brewer.pal(4, "Set1")
names(fill_colors) <- sort(levels(guide_rate$antitag_pos1))

figure_3B <- ggplot(guide_rate, aes(x=antitag_pos1, y=Estimate)) + 
   geom_violin(aes(fill=antitag_pos1), 
               draw_quantiles=c(0.25, 0.75), linetype="dashed") +
   geom_violin(fill="transparent", draw_quantiles=0.5) +
   geom_dotplot(binaxis="y", stackdir="center", binwidth=1) +
   xlab("first position of anti-tag") + ylab("activator-dependent rate\nRFU/min") + 
   theme_classic(base_size=10) + guides(fill="none") + 
   scale_fill_manual(values=fill_colors) + 
   stat_pvalue_manual(data=subset(antitag_pos1_pvalue, group1=="G"),
                      label="p.signif", y.position=c(80, 90, 100), size=2)
   
ggsave(filename=file.path(figure_dir, "figure_3B.pdf"),
       plot=figure_3B,
       device="pdf", width=3, height=2, units="in")
