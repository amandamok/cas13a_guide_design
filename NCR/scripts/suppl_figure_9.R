rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(latex2exp)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data", "supplementary_data")

# generate platemap -------------------------------------------------------

concentrations <- c("no_protein", "no_activator", "250", "125", "63", 
                    "31", "16", "8", "4", "0")
platemap <- data.frame(plate_row=rep(LETTERS[1:4], each=2*length(concentrations)),
                       plate_col=rep(seq(2*length(concentrations)), times=4),
                       conc=rep(rep(concentrations, each=2), times=4))
platemap$plate_well <- with(platemap, paste0(plate_row, plate_col))
platemap$multiplex <- ifelse(platemap$plate_row %in% LETTERS[1:2],
                             "8 crRNAs", "32 crRNAs")
platemap$sample <- with(platemap, paste(multiplex, conc, sep="_"))

# load plate data ---------------------------------------------------------

plate_data_fname <- file.path(data_dir, 
                              "8v32 guide LoD_QC Pairwise Guide Screening (Modified)_20211112_144405.xlsx")
plate_data <- readxl::read_xlsx(plate_data_fname, range="A84:NW144")
plate_data$time <- round(plate_data$`Time [s]`/60)

plate_data <- lapply(platemap$plate_well,
                     function(x) {
                       data.frame(time=plate_data$time,
                                  RFU=plate_data[[x]],
                                  plate_well=x)
                     })
plate_data <- do.call(rbind, plate_data)
plate_data <- plyr::join(plate_data,
                         platemap[, c("plate_well", "multiplex", "conc", "sample")],
                         by="plate_well")

# generate plot -----------------------------------------------------------

figure_S9 <- ggplot(plate_data,
                    aes(x=time, y=RFU, group=plate_well, col=conc)) +
  geom_line() + theme_classic(base_size=8) + xlab("time (min)") + 
  facet_grid(factor(multiplex, levels=unique(platemap$multiplex)) ~ 
               factor(conc, levels=concentrations), scales="free_y") +
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename=file.path(figure_dir, "suppl_figure_9.pdf"),
       plot=figure_S9,
       device="pdf", width=6.5, height=4, units="in")

