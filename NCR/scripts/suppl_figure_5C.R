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

figure_S5C_top <- ggplot(subset(noActivator_rate, 
                                             NCR.id %in% guide_rate$NCR.id), 
                                      aes(x=direct_repeat, y=Estimate)) + 
  geom_dotplot(binaxis="y", stackdir="center", binwidth=1) + 
  xlab("") + ylab("RFU/min") + ggtitle("activator-independent rate") + 
  theme_classic(base_size=8) + theme(axis.text.x=element_blank())

figure_S5C_bottom <- ggplot(guide_rate, aes(x=direct_repeat, y=Estimate)) + 
  geom_dotplot(binaxis="y", stackdir="center", binwidth=1) + 
  xlab("structure of crRNA direct repeat") + ylab("RFU/min") +
  ggtitle("activator-dependent rate") + 
  theme_classic(base_size=8) + theme(axis.text.x=element_text(angle=90, vjust=0))

figure_S5C <- figure_S5C_top / figure_S5C_bottom

ggsave(filename=file.path(figure_dir, "suppl_figure_5C.pdf"),
       plot=figure_S5C,
       device="pdf", width=6.5, height=4, units="in")
