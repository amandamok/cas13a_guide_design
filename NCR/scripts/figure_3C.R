rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data", "supplementary_data")

# load platemap -----------------------------------------------------------

singleLOD_platemap <- openxlsx::read.xlsx(file.path(data_dir, "SingleLOD_2 Platemap.xlsx"),
                                          sheet=1, rows=2:18, cols=2:17)
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

singleLOD_data <- openxlsx::read.xlsx(file.path(data_dir, "SingleLOD_2.xlsx"),
                                      sheet="Sheet2", rows=54:114)
singleLOD_data <- lapply(seq(nrow(singleLOD_platemap)),
                         function(x) {
                           tmp_well <- singleLOD_platemap$plate_well[x]
                           data.frame(time=singleLOD_data$`Time.[s]`,
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
                                     "p < 0.05", "p >= 0.05")

NCR_504_ttest_plot <- ggplot(NCR_504_ttest, 
                             aes(x=time, y=estimate, 
                                 ymin=conf_int_min, ymax=conf_int_max)) + 
  geom_line(col="grey35") + theme_classic(base_size=8) + 
  geom_point(aes(col=pvalue_label), size=1) + geom_hline(yintercept=0) +
  geom_errorbar(aes(col=pvalue_label)) + 
  scale_color_manual(values=c("p >= 0.05"="grey35", "p < 0.05"="red")) + 
  # ylab(bquote("RFU"["100 fM"] - "RFU"["0 fM"])) + 
  xlab("time (min)") + 
  ylab(bquote(atop(Delta*"(reporter cleaved)", "(RFU)"))) +
  theme(legend.position="bottom") + labs(col="")

NCR_504_slope <- lapply(unique(NCR_504_data$time)[-c(1:4)],
                        function(x) {
                          tmp_data <- subset(NCR_504_data, time<=x)
                          # tmp_model <- nlme::lme(RFU_baselineCorr ~ time*conc,
                          #                        random = ~ 1 + time | well, 
                          #                        tmp_data)
                          # tmp_coef <- data.frame(summary(tmp_model)$tTable)
                          tmp_model <- lm(RFU_baselineCorr ~ time*conc, tmp_data)
                          tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                          colnames(tmp_coef)[colnames(tmp_coef)=="Value"] <- "Estimate"
                          colnames(tmp_coef)[colnames(tmp_coef)=="Std..Error"] <- "Std.Error"
                          colnames(tmp_coef)[colnames(tmp_coef)=="Pr...t.."] <- "p.value"
                          tmp_coef$coef <- rownames(tmp_coef)
                          tmp_coef$time <- x
                          return(tmp_coef)
                        })
NCR_504_slope <- do.call(rbind, NCR_504_slope)
NCR_504_slope <- subset(NCR_504_slope, NCR_504_slope$coef == "time:conc100 fM")
NCR_504_slope$pvalue_label <- ifelse(NCR_504_slope$p.value < 0.05,
                                     "p < 0.05", "p >= 0.05")

NCR_504_slope_plot <- ggplot(NCR_504_slope, 
                             aes(x=time, y=Estimate,
                                 ymin=Estimate+qnorm(0.025)*Std.Error,
                                 ymax=Estimate+qnorm(0.975)*Std.Error)) + 
  geom_line(col="grey35") + theme_classic(base_size=8) + 
  geom_point(aes(col=pvalue_label), size=1) + 
  geom_errorbar(aes(col=pvalue_label)) + 
  scale_color_manual(values=c("p >= 0.05"="grey35", "p < 0.05"="red")) + 
  geom_hline(yintercept=0) + xlab("time (min)") + ylab("reaction rate\n(RFU/min)") + 
  theme(legend.position="bottom") + labs(col="")

NCR_504_trace <- ggplot(NCR_504_data, 
                        aes(x=time, y=RFU_baselineCorr, group=well, col=conc)) + 
  geom_line() + geom_point(size=0.25) + theme_classic(base_size=8) +
  scale_color_manual(values=c("100 fM"="red", "0 fM"="grey35")) + 
  labs(color="") + xlab("time (min)") + ylab("RFU") + 
  theme(legend.position="bottom")

figure_3C <- (NCR_504_trace + ggtitle("", subtitle="Raw fluorescence")) + 
  (NCR_504_ttest_plot + ggtitle("", "Endpoint fluorescence")) + 
  (NCR_504_slope_plot + ggtitle("", "Fluorescence rate"))
ggsave(filename=file.path(figure_dir, "figure_3C.pdf"),
       plot=figure_3C,
       device="pdf", width=6.5, height=2.25, units="in")
