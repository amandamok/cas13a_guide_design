rm(list=ls())

library(here)
library(ggplot2)
library(ggpubr)
library(latex2exp)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
suppl_data_dir <- file.path(project_dir, "data", "supplementary_data")

# load platemap -----------------------------------------------------------

platemap <- readxl::read_xlsx(file.path(suppl_data_dir,
                                        "FIND-IT_multiplex_platemap.xlsx"),
                              range="B2:O18")
conc_col_mapping <- seq(ncol(platemap))
names(conc_col_mapping) <- rep(c(125, 63, 31, 16, 8, 4, 0), each=2)

platemap <- data.frame(multiplex=unlist(platemap),
                       plate_row=rep(LETTERS[seq(nrow(platemap))],
                                     times=ncol(platemap)),
                       plate_col=rep(seq(ncol(platemap)),
                                     each=nrow(platemap)))
platemap <- subset(platemap, !is.na(platemap$multiplex))
platemap$multiplex <- paste(platemap$multiplex, "crRNAs")
platemap$plate_well <- with(platemap, paste0(plate_row, plate_col))
platemap$conc <- names(conc_col_mapping)[platemap$plate_col]
platemap$conc <- factor(platemap$conc, levels=unique(names(conc_col_mapping)))

# load plate reader data --------------------------------------------------

plate_data <- readxl::read_xlsx(file.path(suppl_data_dir,
                                          "FIND-IT_multiplex.xlsx"),
                                range="A82:NW142")
plate_data$time <- plate_data$`Time [s]`/60

plate_data <- lapply(seq(nrow(platemap)),
                     function(x) {
                       tmp_time <- plate_data$time
                       f6_timepoint <- which.min(abs(tmp_time-6))
                       tmp_RFU <- plate_data[[platemap$plate_well[x]]]
                       f6_RFU <- tmp_RFU[f6_timepoint]
                       norm_RFU <- tmp_RFU / f6_RFU
                       return(data.frame(time=tmp_time,
                                         RFU=norm_RFU,
                                         plate_well=platemap$plate_well[x],
                                         multiplex=platemap$multiplex[x],
                                         conc=platemap$conc[x]))
                     })
plate_data <- do.call(rbind, plate_data)

trace_plots <- ggplot(plate_data, 
                      aes(x=time, y=RFU, group=plate_well,
                          col=factor(conc, levels=unique(platemap$conc)))) + 
  geom_line() + facet_grid(~multiplex) + theme_classic() + labs(col="conc [cp/uL]")

# test endpoint fluorescence ----------------------------------------------

endpoint_data <- subset(plate_data, time==max(plate_data$time))

endpoint_pvalues <- lapply(unique(platemap$multiplex),
                           function(x) {
                             within(rstatix::t_test(subset(endpoint_data,
                                                           multiplex==x),
                                                    RFU ~ conc, ref.group="0",
                                                    paired=F, var.equal=F,
                                                    alternative="less",
                                                    p.adjust.method="none"),
                                    multiplex <- x)
                           })
endpoint_pvalues <- do.call(rbind, endpoint_pvalues)
endpoint_pvalues$y.position <- rep(c(7.3, 5.5, 5.2, 4.9, 4.6, 4.3), times=2)

# generate plot -----------------------------------------------------------

endpoint_means <- within(aggregate(RFU ~ conc + multiplex, 
                                   endpoint_data, FUN=mean),
                         std_dev <- aggregate(RFU ~ conc + multiplex, 
                                              endpoint_data, FUN=sd)$RFU)
endpoint_means$pvalue <- endpoint_pvalues$p[prodlim::row.match(endpoint_means[, c("conc", "multiplex")],
                                                               endpoint_pvalues[, c("group2", "multiplex")])]
endpoint_means$pvalue_label <- ifelse(endpoint_means$pvalue < 0.05,
                                      "p < 0.05", "ns")
endpoint_means$pvalue_label[is.na(endpoint_means$pvalue_label)] <- "ns"

ggplot(endpoint_means, aes(x=conc, y=RFU)) + 
  geom_col(aes(fill=pvalue_label)) + theme_classic(base_size=10) + 
  geom_errorbar(aes(ymin=RFU + qnorm(0.025)*std_dev,
                    ymax=RFU + qnorm(0.975)*std_dev), width=0.5) + 
  facet_grid(~factor(multiplex, levels=c("8 crRNAs", "32 crRNAs"))) + 
  geom_point(data=endpoint_data, shape=1) + 
  stat_pvalue_manual(endpoint_pvalues) + labs(fill="") + 
  scale_fill_manual(values=c("p < 0.05"="red", "ns"="grey35")) + 
  xlab(TeX("SARS-CoV-2 RNA (copies/$\\mu$L)")) + 
  ylab(TeX("Fluorescence ($F/F_6$)")) + 
  theme(panel.spacing.x=unit(3, "lines"))
