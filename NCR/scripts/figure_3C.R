rm(list=ls())

library(here)
library(ggplot2)
library(ggpubr)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))

# calculate differences ---------------------------------------------------

guide_rate$antitag_num_complementary <- sapply(guide_rate$antitag,
                                               function(x) {
                                                 tmp_antitag <- strsplit(x, split="")[[1]]
                                                 complementary_antitag <- c("G", "U", "U", "U")
                                                 sum(sapply(seq_along(tmp_antitag),
                                                            function(y) {
                                                              tmp_antitag[y] == complementary_antitag[y]
                                                            }))
                                               })
guide_rate$antitag_num_complementary <- factor(guide_rate$antitag_num_complementary)
levels(guide_rate$antitag_num_complementary) <- sapply(levels(guide_rate$antitag_num_complementary),
                                                       function(x) {
                                                         num_guides <- sum(guide_rate$antitag_num_complementary == x)
                                                         tmp_label <- paste0(x, "\n", "(n=", num_guides, ")")
                                                         return(tmp_label)
                                                       })
num_complementary_pvalue <- compare_means(Estimate ~ antitag_num_complementary,
                                          guide_rate, method="wilcox.test", 
                                          p.adjust.method="none")

# generate plots ----------------------------------------------------------

figure_3C <- ggplot(guide_rate, aes(x=antitag_num_complementary, y=Estimate)) + 
  geom_violin(aes(fill=antitag_num_complementary), 
              draw_quantiles=c(0.25, 0.75), linetype="dashed") +
  geom_violin(fill="transparent", draw_quantiles=0.5) +
  # geom_dotplot(binaxis="y", stackdir="center", binwidth=1) + 
  xlab("# complementary positions") + ylab("activator-dependent rate\nRFU/min") + 
  theme_classic(base_size=10) + guides(fill="none") + 
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Blues")) + 
  stat_pvalue_manual(data=subset(num_complementary_pvalue, group1=="0\n(n=47)"),
                     label="p.signif", y.position=c(80, 90, 100), size=2)

ggsave(filename=file.path(figure_dir, "figure_3C.pdf"),
       plot=figure_3C,
       device="pdf", width=3, height=2, units="in")
