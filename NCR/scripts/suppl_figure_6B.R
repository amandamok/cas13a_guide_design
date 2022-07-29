themerm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load viral genome annotations -------------------------------------------

viral_genome <- read.csv(file.path(data_dir, "viral_genome_features.csv"))

# load computed slopes ----------------------------------------------------

guide_rate <- read.csv(file.path(data_dir, "guide_rate.csv"))

# parse seed regions (spacer pos 9-13) ------------------------------------

seed_start <- 9
seed_stop <- 13

target_start <- c(20:1)[seed_stop]
target_stop <- c(20:1)[seed_start]

guide_structure_seed <- sapply(seq(nrow(guide_rate)),
                               function(x) {
                                 region_start <- guide_rate$start[x] + target_start - 1
                                 region_end <- guide_rate$start[x] + target_stop - 1
                                 region_features <- viral_genome[match(region_start:region_end,
                                                                       viral_genome$pos),]
                                 region_annotations <- c(manfredonia = mean(region_features$manfredonia == "lowShannon/highSHAPE"),
                                                         lan = mean(region_features$lan == "unstructured"),
                                                         sun_invivo = mean(region_features$sun_invivo),
                                                         sun_invitro = mean(region_features$sun_invitro),
                                                         huston = mean(region_features$huston == 0))
                                 return(region_annotations)
                               })
guide_structure_seed <- data.frame(t(guide_structure_seed))
colnames(guide_structure_seed) <- paste0(colnames(guide_structure_seed), "_seed")
guide_rate <- cbind(guide_rate, guide_structure_seed)

# generate plots ----------------------------------------------------------

lan_target <- ggplot(guide_rate, aes(x=lan, y=Estimate)) + 
  theme_classic(base_size=8) + geom_point(shape=16, col="grey35") + 
  geom_smooth(method=lm, formula=y~x) +
  xlab("% unstructured") + ylab("activator-dependent rate\nRFU/min") + 
  ggtitle("Lan (2020)", subtitle="20nt target region")

lan_seed <- ggplot(guide_rate, aes(x=lan_seed, y=Estimate)) + 
  theme_classic(base_size=8) + geom_point(shape=16, col="grey35") + 
  geom_smooth(method=lm, formula=y~x) +
  xlab("% unstructured") + ylab("activator-dependent rate\nRFU/min") + 
  ggtitle("", subtitle="target seed region")

huston_target <- ggplot(guide_rate, aes(x=huston, y=Estimate)) + 
  theme_classic(base_size=8) + geom_point(shape=16, col="grey35") + 
  geom_smooth(method=lm, formula=y~x) +
  xlab("% double-stranded") + ylab("") +
  ggtitle("Huston (2021)", subtitle="20nt target region")

huston_seed <- ggplot(guide_rate, aes(x=huston_seed, y=Estimate)) + 
  theme_classic(base_size=8) + geom_point(shape=16, col="grey35") + 
  geom_smooth(method=lm, formula=y~x) +
  xlab("% double-stranded") + ylab("") +
  ggtitle("", subtitle="target seed region")

sun_invivo_target <- ggplot(guide_rate, aes(x=sun_invivo, y=Estimate)) + 
  theme_classic(base_size=8) + geom_point(shape=16, col="grey35") + 
  geom_smooth(method=lm, formula=y~x) +
  xlab("mean icSHAPE score") + ylab("") +
  ggtitle("Sun (2021)", subtitle="20nt target region")

sun_invivo_seed <- ggplot(guide_rate, aes(x=sun_invivo_seed, y=Estimate)) + 
  theme_classic(base_size=8) + geom_point(shape=16, col="grey35") + 
  geom_smooth(method=lm, formula=y~x) +
  xlab("mean icSHAPE score") + ylab("") +
  ggtitle("", subtitle="target seed region")

suppl_figure_6B <- lan_target + huston_target + sun_invivo_target + 
  lan_seed + huston_seed + sun_invivo_seed

ggsave(file=file.path(figure_dir, "suppl_figure_6B.pdf"),
       plot=suppl_figure_6B,
       device="pdf", width=6.5, height=4, units="in")
