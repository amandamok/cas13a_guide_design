rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load computed slopes ----------------------------------------------------

guide_rate <- read.csv(file.path(data_dir, "guide_rate.csv"))
guide_rate$direct_repeat <- relevel(as.factor(guide_rate$direct_repeat),
                                    ref=".......((((.........))))........")

noActivator_rate <- read.csv(file.path(data_dir, "noActivator_rate.csv"))
noActivator_rate$direct_repeat <- factor(noActivator_rate$direct_repeat,
                                         levels=levels(guide_rate$direct_repeat))

# generate plots ----------------------------------------------------------

suppl_figure_4B_noActivator <- ggplot(subset(noActivator_rate, 
                                             NCR.id %in% guide_rate$NCR.id), 
                                      aes(x=direct_repeat, y=Estimate)) + 
  geom_dotplot(binaxis="y", stackdir="center", binwidth=1) + 
  xlab("") + ylab("RFU/min") + ggtitle("activator-independent rate") + 
  theme_classic(base_size=10) + theme(axis.text.x=element_blank())

suppl_figure_4B_100fM <- ggplot(guide_rate, aes(x=direct_repeat, y=Estimate)) + 
  geom_dotplot(binaxis="y", stackdir="center", binwidth=1) + 
  xlab("structure of crRNA direct repeat") + ylab("RFU/min") +
  ggtitle("activator-dependent rate") + 
  theme_classic(base_size=10) + theme(axis.text.x=element_text(angle=90, vjust=0))

suppl_figure_4B <- suppl_figure_4B_noActivator / suppl_figure_4B_100fM

ggsave(filename=file.path(figure_dir, "suppl_figure_4B.pdf"),
       plot=suppl_figure_4B,
       device="pdf", width=6.5, height=5, units="in")
