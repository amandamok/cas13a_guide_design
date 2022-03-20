rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load computed slopes ----------------------------------------------------

guide_rate <- read.csv(file.path(data_dir, "guide_rate.csv"))

# generate plots ----------------------------------------------------------

lan_huston <- ggplot(guide_rate, aes(x=lan, y=huston)) + 
  geom_point(col="grey35") + geom_smooth(method="lm", formula=y~x) + 
  theme_classic(base_size=10) + xlab("Lan (2020)") + ylab("Huston (2021)") + 
  ggtitle(bquote(rho == .(with(guide_rate, 
                               signif(cor(lan, huston, method="spearman"), 3)))))

sun_invivo_huston <- ggplot(guide_rate, aes(x=sun_invivo, y=huston)) + 
  geom_point(col="grey35") + geom_smooth(method="lm", formula=y~x) + 
  theme_classic(base_size=10) + xlab("Sun (2021): in vivo") + ylab("Huston (2021)") + 
  ggtitle(bquote(rho == .(with(guide_rate, 
                               signif(cor(sun_invivo, huston, method="spearman"), 3)))))

lan_sun_invivo <- ggplot(guide_rate, aes(x=lan, y=sun_invivo)) + 
  geom_point(col="grey35") + geom_smooth(method="lm", formula=y~x) + 
  theme_classic(base_size=10) + xlab("Lan (2020)") + ylab("Sun (2021): in vivo") + 
  ggtitle(bquote(rho == .(with(guide_rate, 
                               signif(cor(lan, sun_invivo, method="spearman"), 3)))))

suppl_figure_5A <- lan_huston + sun_invivo_huston + lan_sun_invivo

ggsave(file=file.path(figure_dir, "suppl_figure_5A.pdf"),
       plot=suppl_figure_5A,
       device="pdf", width=6.5, height=2, units="in")
