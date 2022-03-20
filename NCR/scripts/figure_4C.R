rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data", "supplementary_data")

# load platemap -----------------------------------------------------------

singleLOD_platemap <- xlsx::read.xlsx(file.path(data_dir, "SingleLOD_2 Platemap.xlsx"),
                                      sheetIndex=1, startRow=2, endRow=18, colIndex=2:17)
singleLOD_platemap <- data.frame(plate_row=rep(LETTERS[1:nrow(singleLOD_platemap)],
                                               times=ncol(singleLOD_platemap)),
                                 plate_col=rep(seq(ncol(singleLOD_platemap)),
                                               each=nrow(singleLOD_platemap)),
                                 guide_id=unlist(singleLOD_platemap),
                                 conc=rep(c("No protein", "No activator",
                                            "1 pM", "300 fM", "100 fM",
                                            "30 fM", "10 fM", "3 fM"),
                                          each=2*nrow(singleLOD_platemap)))
singleLOD_platemap$plate_well <- with(singleLOD_platemap, paste0(plate_row, plate_col))
singleLOD_platemap <- subset(singleLOD_platemap, guide_id != "NA")

plate_guides <- unique(singleLOD_platemap$guide_id)
plate_conc <- unique(singleLOD_platemap$conc)

# load data ---------------------------------------------------------------

singleLOD_data <- xlsx::read.xlsx(file.path(data_dir, "SingleLOD_2.xlsx"),
                                  sheetName="Sheet2", startRow=54, endRow=114)
singleLOD_data <- lapply(seq(nrow(singleLOD_platemap)),
                         function(x) {
                           tmp_well <- singleLOD_platemap$plate_well[x]
                           data.frame(time=singleLOD_data$Time..s.,
                                      RFU=as.numeric(singleLOD_data[[tmp_well]]),
                                      well=tmp_well,
                                      guide_id=singleLOD_platemap$guide_id[x],
                                      conc=singleLOD_platemap$conc[x])
                         })
singleLOD_data <- do.call(rbind, singleLOD_data)
singleLOD_data$time <- round(singleLOD_data$time / 60)

# baseline correction -----------------------------------------------------

baseline_RFU <- subset(singleLOD_data, time==10)
singleLOD_data$baseline_RFU <- baseline_RFU$RFU[match(singleLOD_data$well, 
                                                      baseline_RFU$well)]
singleLOD_data$RFU_baselineCorr <- with(singleLOD_data, RFU-baseline_RFU)

# compute slopes ----------------------------------------------------------

all_slopes <- lapply(plate_guides,
                     function(x) {
                       tmp_slopes <- lapply(plate_conc,
                                            function(y) {
                                              tmp_data <- subset(singleLOD_data, 
                                                                 guide_id==x &
                                                                   conc==y & 
                                                                   time >= 10)
                                              tmp_model <- lm(RFU ~ time, tmp_data)
                                              tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                              tmp_coef$coef <- rownames(tmp_coef)
                                              tmp_coef$guide_id <- x
                                              tmp_coef$conc <- y
                                              return(tmp_coef)
                                            })
                       tmp_slopes <- do.call(rbind, tmp_slopes)
                       return(tmp_slopes)
                     })
all_slopes <- do.call(rbind, all_slopes)
all_slopes <- subset(all_slopes, coef=="time")
all_slopes$conc <- factor(all_slopes$conc, levels=unique(all_slopes$conc))

# compare slopes to "no activator" ----------------------------------------

LOD_pvalues <- lapply(plate_guides,
                      function(x) {
                        tmp_pvalues <- lapply(plate_conc[-c(1,2)],
                                              function(y) {
                                                tmp_data <- subset(singleLOD_data,
                                                                   guide_id==x & 
                                                                     time >= 10 &
                                                                     conc %in% c(y, "No activator"))
                                                tmp_data$conc <- factor(tmp_data$conc,
                                                                        levels=c("No activator", y))
                                                tmp_model <- lm(RFU~time*conc, tmp_data)
                                                tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                tmp_coef$coef <- rownames(tmp_coef)
                                                tmp_coef$guide_id <- x
                                                tmp_coef$conc <- y
                                                return(tmp_coef)
                                              })
                        tmp_pvalues <- do.call(rbind, tmp_pvalues)
                        return(tmp_pvalues)
                      })
LOD_pvalues <- do.call(rbind, LOD_pvalues)
LOD_pvalues <- subset(LOD_pvalues, grepl("time:", LOD_pvalues$coef))

# annotate all_slopes with LOD pvalues
all_slopes$LOD_pvalue <- LOD_pvalues$Pr...t..[prodlim::row.match(all_slopes[, c("guide_id", "conc")],
                                                                 LOD_pvalues[, c("guide_id", "conc")])]
all_slopes$LOD_pvalue_label <- ifelse(all_slopes$LOD_pvalue < 0.05, "p<0.05", "p≥0.05")

# generate LOD plots ------------------------------------------------------

LOD_plots <- ggplot(all_slopes, aes(x=conc, y=Estimate,
                                    ymin=Estimate+qnorm(0.025)*Std..Error,
                                    ymax=Estimate+qnorm(0.975)*Std..Error,
                                    fill=LOD_pvalue_label)) + 
  geom_col() + geom_errorbar() + theme_classic() + 
  facet_wrap(vars(guide_id), scales="free_y", nrow=2) + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
  xlab("") + ylab("RFU/min") + labs(fill=" ") + 
  scale_fill_manual(values=c("p<0.05" = "red", "p≥0.05" = "grey35"))

# endpoint vs. slope: NCR_504 @ 100fM -------------------------------------

NCR_504_data <- subset(singleLOD_data, guide_id == 504 & 
                         time >= 10 & conc %in% c("100 fM", "No activator"))
NCR_504_data$conc[NCR_504_data$conc=="No activator"] <- "0 fM"
NCR_504_data$conc <- factor(NCR_504_data$conc, levels=c("0 fM", "100 fM"))

NCR_504_ttest <- lapply(unique(NCR_504_data$time)[-1],
                        function(x) {
                          tmp_data <- subset(NCR_504_data, time==x)
                          tmp_data$conc <- relevel(tmp_data$conc, ref="100 fM")
                          tmp_test <- t.test(RFU_baselineCorr ~ conc, tmp_data)
                          data.frame(t_statistic=tmp_test$statistic,
                                     df=tmp_test$parameter,
                                     pvalue=tmp_test$p.value,
                                     estimate=tmp_test$estimate[1] - tmp_test$estimate[2],
                                     std_err=tmp_test$stderr,
                                     conf_int_min=tmp_test$conf.int[1],
                                     conf_int_max=tmp_test$conf.int[2],
                                     time=x)
                        })
NCR_504_ttest <- do.call(rbind, NCR_504_ttest)
NCR_504_ttest$pvalue_label <- ifelse(NCR_504_ttest$pvalue < 0.05,
                                     "p<0.05", "p>=0.05")

NCR_504_ttest_plot <- ggplot(NCR_504_ttest, 
                             aes(x=time, y=estimate, 
                                 ymin=conf_int_min, ymax=conf_int_max)) + 
  geom_point(aes(col=pvalue_label), size=0.5) + 
  geom_errorbar(aes(col=pvalue_label)) + 
  geom_line(col="grey35") + theme_classic(base_size=10) + 
  scale_color_manual(values=c("p>=0.05"="grey35", "p<0.05"="red")) + 
  geom_hline(yintercept=0) + ylab(bquote("RFU"["100 fM"] - "RFU"["0 fM"])) + 
  theme(legend.position="top") + labs(col="")

NCR_504_slope <- lapply(unique(NCR_504_data$time)[-1],
                        function(x) {
                          tmp_data <- subset(NCR_504_data, time<=x)
                          tmp_model <- lm(RFU_baselineCorr ~ time*conc,
                                          tmp_data)
                          tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                          tmp_coef$coef <- rownames(tmp_coef)
                          tmp_coef$time <- x
                          return(tmp_coef)
                        })
NCR_504_slope <- do.call(rbind, NCR_504_slope)
NCR_504_slope <- subset(NCR_504_slope, NCR_504_slope$coef == "time:conc100 fM")
NCR_504_slope$pvalue_label <- ifelse(NCR_504_slope$Pr...t.. < 0.05,
                                     "p<0.05", "p>=0.05")

NCR_504_slope_plot <- ggplot(NCR_504_slope, 
                             aes(x=time, y=Estimate,
                                 ymin=Estimate+qnorm(0.025)*Std..Error,
                                 ymax=Estimate+qnorm(0.975)*Std..Error)) + 
  geom_point(aes(col=pvalue_label), size=0.5) + 
  geom_errorbar(aes(col=pvalue_label)) + 
  geom_line(col="grey35") + theme_classic(base_size=10) + 
  scale_color_manual(values=c("p>=0.05"="grey35", "p<0.05"="red")) + 
  geom_hline(yintercept=0) + ylab("RFU/min") + 
  theme(legend.position="top") + labs(col="")

NCR_504_trace <- ggplot(NCR_504_data, 
       aes(x=time, y=RFU_baselineCorr, group=well, col=conc)) + 
  geom_line() + geom_point(size=0.25) + theme_classic(base_size=10) +
  scale_color_manual(values=c("100 fM"="red", "0 fM"="grey35")) + 
  labs(color="") + xlab("min") + ylab("RFU") + 
  theme(legend.position="top")

figure_4C <- NCR_504_trace + NCR_504_ttest_plot + NCR_504_slope_plot
ggsave(filename=file.path(figure_dir, "figure_4C.pdf"),
       plot=figure_4C,
       device="pdf", width=6.5, height=2.5, units="in")
