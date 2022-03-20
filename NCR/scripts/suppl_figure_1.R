rm(list=ls())

library(here)
library(ggplot2)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load computed slopes ----------------------------------------------------

guide_rate <- read.csv(file.path(data_dir, "guide_rate.csv"))
gblock_rate <- read.csv(file.path(data_dir, "gblock_rate.csv"))

# merge datasets ----------------------------------------------------------

guide_rate$gblock_rate <- gblock_rate$Estimate[match(guide_rate$guide_id,
                                                     gblock_rate$guide_id)]
guide_rate$gblock_stderr <- gblock_rate$Std.Error[match(guide_rate$guide_id,
                                                        gblock_rate$guide_id)]

# remove outlier
guide_rate <- subset(guide_rate, guide_id != "1352")

# remove 28mers

guide_rate <- subset(guide_rate, nchar(spacer) == 20)

# generate plot -----------------------------------------------------------

vRNA_gBlock_corr <- with(guide_rate, cor(Estimate, gblock_rate))

suppl_figure_1 <- ggplot(guide_rate, 
                         aes(x=Estimate, y=gblock_rate,
                             xmin=Estimate+qnorm(0.025)*Std.Error,
                             xmax=Estimate+qnorm(0.975)*Std.Error,
                             ymin=gblock_rate+qnorm(0.025)*gblock_stderr,
                             ymax=gblock_rate+qnorm(0.975)*gblock_stderr)) + 
  geom_errorbar(col="grey35", alpha=0.5) + geom_errorbarh(col="grey35", alpha=0.5) + 
  geom_point(col="grey35") + theme_classic(base_size=10) + coord_fixed(ratio=1) +
  xlab("vRNA rate (RFU/min)") + ylab("gBlock rate (RFU/min)") + 
  ggtitle(bquote(rho == .(signif(vRNA_gBlock_corr, 3))))

ggsave(filename=file.path(figure_dir, "suppl_figure_1.pdf"),
       plot=suppl_figure_1,
       device="pdf", width=5, units="in")
