rm(list=ls())

library(here)
library(ggplot2)
library(ggpubr)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))

# calculate differences ---------------------------------------------------

guide_rate$antitag_group <- guide_rate$antitag_label
guide_rate$antitag_group[!grepl("G", guide_rate$antitag_label)] <- "non-G"
guide_rate$antitag_group <- factor(guide_rate$antitag_group,
                                   levels=c("non-G", "G", "GU", "GUU"))
levels(guide_rate$antitag_group) <- sapply(levels(guide_rate$antitag_group),
                                           function(x) {
                                             num_guides <- sum(guide_rate$antitag_group == x)
                                             tmp_label <- paste0(x, "\n", "(n=", num_guides, ")")
                                             return(tmp_label)
                                           })
antitag_group_pvalue <- compare_means(Estimate ~ antitag_group,
                                      guide_rate, method="wilcox.test",
                                      p.adjust.method="none")

# generate plots ----------------------------------------------------------

figure_3D <- ggplot(guide_rate, aes(x=antitag_group, y=Estimate)) + 
  geom_violin(aes(fill=antitag_group),
              draw_quantiles=c(0.25, 0.75), linetype="dashed") + 
  geom_violin(fill="transparent", draw_quantiles=0.5) + 
  # geom_dotplot(binaxis="y", stackdir="center", binwidth=1) + 
  xlab("anti-tag") + ylab("activator-dependent rate\nRFU/min") + 
  theme_classic(base_size=10) + guides(fill="none") + 
  stat_pvalue_manual(data=subset(antitag_group_pvalue, group1=="non-G\n(n=164)"),
                     label="p.signif", y.position=c(85, 95, 105), size=2) + 
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "BuGn"))

ggsave(filename=file.path(figure_dir, "figure_3D.pdf"),
       plot=figure_3D,
       device="pdf", width=3, height=2, units="in")
