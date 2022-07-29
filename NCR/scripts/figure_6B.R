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

# calculate slopes --------------------------------------------------------

lod_slopes <- lapply(unique(platemap$sample),
                     function(x) {
                       tmp_data <- subset(plate_data,
                                          time >= 10 & sample == x)
                       tmp_fit <- lm(RFU ~ time, tmp_data)
                       tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                       tmp_coef$term <- rownames(tmp_coef)
                       tmp_coef$sample <- x
                       return(tmp_coef)
                     })
lod_slopes <- do.call(rbind, lod_slopes)
lod_slopes <- subset(lod_slopes, term=="time")
lod_slopes$multiplex <- sub("_.*$", "", lod_slopes$sample)
lod_slopes$conc <- sub("^.*crRNAs_", "", lod_slopes$sample)

# calculate pvalues -------------------------------------------------------

lod_pvalues <- lapply(unique(platemap$multiplex),
                      function(tmp_multiplex) {
                        tmp_p <- lapply(concentrations[!(concentrations %in% c("no_protein", "no_activator", "0"))],
                                        function(tmp_conc) {
                                          tmp_data <- subset(plate_data,
                                                             time >= 10 & 
                                                               multiplex == tmp_multiplex & 
                                                               conc %in% c("0", tmp_conc))
                                          tmp_fit <- lm(RFU ~ time*conc, tmp_data)
                                          tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                          tmp_coef$term <- rownames(tmp_coef)
                                          tmp_coef$multiplex <- tmp_multiplex
                                          tmp_coef$conc <- tmp_conc
                                          return(tmp_coef)
                                        })
                        tmp_p <- do.call(rbind, tmp_p)
                        return(tmp_p)
                      })
lod_pvalues <- do.call(rbind, lod_pvalues)
lod_pvalues <- subset(lod_pvalues, grepl("time:", lod_pvalues$term))

# for stat_pvalue_manual: group1, group2, p, multiplex, y.position, label
lod_pvalues$.y. <- "RFU"
lod_pvalues$group1 <- "0"
lod_pvalues$group2 <- rep(concentrations[!(concentrations %in% c("no_protein", "no_activator", "0"))],
                          times=2)
lod_pvalues$p.label <- cut(lod_pvalues$Pr...t..,
                           breaks=c(0, 1e-5, 1e-4, 1e-3, 0.05, 1),
                           labels=c("****", "***", "**", "*", "ns"))
lod_pvalues$y.position <- c(seq(from=20, to=7.5, length.out=7),
                            seq(from=25, to=7.5, length.out=7))

lod_slopes$LOD_p <- lod_pvalues$Pr...t..[prodlim::row.match(lod_slopes[, c("multiplex", "conc")],
                                                            lod_pvalues[, c("multiplex", "conc")])]
lod_slopes$LOD_p_label <- ifelse(lod_slopes$LOD_p < 0.05,
                                 "p < 0.05", "ns")
lod_slopes$LOD_p_label[is.na(lod_slopes$LOD_p_label)] <- "ns"

# generate plot -----------------------------------------------------------

figure_6B <- ggplot(subset(lod_slopes, !grepl("no_", conc)),
                    aes(x=factor(conc, levels=concentrations[-c(1,2)]), y=Estimate, 
                        ymin=Estimate + qnorm(0.025)*Std..Error,
                        ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col(aes(fill=LOD_p_label)) + geom_errorbar(width=0.5) + 
  facet_grid(~factor(multiplex, levels=c("8 crRNAs", "32 crRNAs"))) + 
  theme_classic(base_size=8) + stat_pvalue_manual(lod_pvalues, label="p.label") + 
  scale_fill_manual(values=c("p < 0.05"="red", "ns"="grey35")) + 
  xlab(TeX("SARS-CoV-2 RNA (copies/$\\mu$L)")) + 
  ylab("RFU/min") + 
  theme(panel.spacing.x=unit(3, "lines")) + labs(fill="")

ggsave(filename=file.path(figure_dir, "figure_6B.pdf"),
       plot=figure_6B,
       device="pdf", width=6.5, height=2, units="in")

