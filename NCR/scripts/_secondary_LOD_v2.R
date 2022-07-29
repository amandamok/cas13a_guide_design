rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(latex2exp)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data", "supplementary_data")

# generate plate map ------------------------------------------------------

data_fname <- file.path(data_dir, "EM_Cas13Lod_Rep3_x8_06152022_18966Zh_80of32gain.xlsx")

concentrations <- c("no_protein", "no_activator",
                    "125", "63", "31", "16", "8", "4", "2", "0")
platemap <- data.frame(plate_row=rep(LETTERS[5:12], each=20),
                       plate_col=rep(c(3:22), times=8),
                       conc=rep(rep(concentrations, each=2), times=8))
platemap$multiplex <- ifelse(platemap$plate_row %in% LETTERS[5:8], "8g", "32g")
platemap$plate_well <- with(platemap, paste(plate_row, plate_col, sep=""))
platemap$sample <- with(platemap, paste(multiplex, conc, sep="_"))

# load plate reader data --------------------------------------------------

plate_data <- readxl::read_xlsx(data_fname, range="A76:NW135")
plate_data <- reshape2::melt(plate_data, id.vars=1:3, variable.name="plate_well")
plate_data$time <- plate_data$`Time [s]`/60
plate_data$RFU <- as.numeric(plate_data$value)
plate_data <- subset(plate_data, plate_well %in% platemap$plate_well)
plate_data <- droplevels(plate_data)
plate_data <- plyr::join(plate_data, platemap, by="plate_well")

# calculate slope by plate well -------------------------------------------

individual_slopes <- lapply(platemap$plate_well,
                            function(x) {
                              tmp_data <- subset(plate_data, 
                                                 plate_well == x & time >= 10)
                              tmp_fit <- lm(RFU ~ time, tmp_data)
                              tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                              tmp_coef$term <- rownames(tmp_coef)
                              tmp_coef$plate_well <- x
                              return(tmp_coef)
                            })
individual_slopes <- do.call(rbind, individual_slopes)
individual_slopes <- subset(individual_slopes, term=="time")
individual_slopes <- plyr::join(individual_slopes, platemap, by="plate_well")

# select middle 4 replicates ----------------------------------------------

which_wells <- lapply(split(individual_slopes, individual_slopes$sample),
                      function(x) {
                        x <- x[order(x$Estimate),]
                        return(x$plate_well[3:6])
                      })
which_wells <- unlist(which_wells)

plate_data_filtered <- subset(plate_data, plate_well %in% which_wells)

# calculate LOD -----------------------------------------------------------

lod_slopes <- lapply(unique(plate_data_filtered$sample),
                     function(x) {
                       tmp_data <- subset(plate_data_filtered,
                                          sample==x & time >= 10)
                       tmp_fit <- lm(RFU ~ time, tmp_data)
                       tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                       tmp_coef$term <- rownames(tmp_coef)
                       tmp_coef$sample <- x
                       return(tmp_coef)
                     })
lod_slopes <- do.call(rbind, lod_slopes)
lod_slopes <- subset(lod_slopes, term=="time")
lod_slopes$multiplex <- factor(sub("_.*$", "", lod_slopes$sample),
                               levels=c("8g", "32g"))
lod_slopes$conc <- factor(sub("^.*g_", "", lod_slopes$sample),
                          levels=concentrations)

lod_pvalues <- lapply(concentrations[-c(1,2)],
                      function(tmp_conc) {
                        tmp_p <- lapply(c("8g", "32g"),
                                        function(tmp_multiplex) {
                                          tmp_data <- subset(plate_data_filtered,
                                                             time >= 10 &
                                                               multiplex == tmp_multiplex & 
                                                               conc %in% c("no_activator", tmp_conc))
                                          tmp_data$conc <- factor(tmp_data$conc,
                                                                  levels=c("no_activator", tmp_conc))
                                          tmp_fit <- lm(RFU ~ time*conc, tmp_data)
                                          tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                          tmp_coef$term <- rownames(tmp_coef)
                                          tmp_coef$conc <- tmp_conc
                                          tmp_coef$multiplex <- tmp_multiplex
                                          return(tmp_coef)
                                        })
                        tmp_p <- do.call(rbind, tmp_p)
                        return(tmp_p)
                      })
lod_pvalues <- do.call(rbind, lod_pvalues)
lod_pvalues <- subset(lod_pvalues, grepl("time:", lod_pvalues$term))

lod_slopes$LOD_p <- lod_pvalues$Pr...t..[prodlim::row.match(lod_slopes[, c("conc", "multiplex")],
                                                            lod_pvalues[, c("conc", "multiplex")])]
lod_slopes$p_label <- ifelse(lod_slopes$LOD_p < 0.05,
                             "p < 0.05", "ns")

ggplot(lod_slopes, aes(x=conc, y=Estimate, fill=p_label,
                       ymin=Estimate + qnorm(0.025)*Std..Error,
                       ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col() + geom_errorbar() + 
  facet_grid(~multiplex) + theme_classic() + 
  scale_fill_manual(values=c("p < 0.05"="red", "ns"="darkgrey")) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
  xlab("") + ylab("RFU/min") + 
  ggtitle("middle 4 replicates")
