rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data")

time_start <- 0
time_stop <- 120

# load plate data ---------------------------------------------------------

load(file.path(data_dir, "vRNA_plate_data.Rda"))

# good guide: NCR_527 -----------------------------------------------------

data_527 <- subset(plate_data, guide_id=="527")
data_527$activator <- relevel(data_527$activator, ref="noActivator")
# compute fit
fit_527 <- lm(RFU ~ time*activator, subset(data_527, time >= 10))
# predict trajectory
data_predict_527 <- data.frame(time=rep(seq(time_start, time_stop, 1), times=2),
                               activator=rep(c("noActivator", "100fM"), 
                                             each=length(seq(time_start, time_stop, 1))))
data_predict_527$RFU <- predict(fit_527, data_predict_527)
# pull regression coefficients
coef_527 <- data.frame(summary(fit_527)$coefficients)
# generate plot text
max_527 <- max(data_527$RFU, na.rm=T)
width_527 <- max(data_527$RFU, na.rm=T) - min(data_527$RFU, na.rm=T)
rate_noActivator_527 <- coef_527$Estimate[rownames(coef_527)=="time"]
rate_100fM_527 <- coef_527$Estimate[rownames(coef_527)=="time"] + 
  coef_527$Estimate[rownames(coef_527)=="time:activator100fM"]
pvalue_527 <- coef_527$Pr...t..[rownames(coef_527)=="time:activator100fM"]
plot_text_527 <- data.frame(x=10,
                            y=max_527 - c(1, 2.5, 4)*width_527/20,
                            text=c(paste("RNP-only rate:", signif(rate_noActivator_527, 3)),
                                   paste("100 fM rate:", signif(rate_100fM_527, 3)),
                                   paste("p =", signif(pvalue_527, 3))))
# generate plot
guide_plot_527 <- ggplot(data_predict_527, aes(x=time, y=RFU, col=activator)) +
  geom_line(size=2) + theme_classic(base_size=10) + 
  geom_point(data=data_527, shape=16, alpha=0.75,
             aes(x=time, y=RFU, col=activator, group=well_384)) +
  scale_color_manual(values=c("noActivator"="grey35", "100fM"="red")) +
  ggtitle("target region: 7720-7739",
          subtitle=paste(plot_text_527$text, collapse="\n")) + 
  xlab("time (min)") + ylab("RFU") + guides(col="none")

# bad guide: NCR_505 ------------------------------------------------------

data_505 <- subset(plate_data, guide_id=="505")
data_505$activator <- relevel(data_505$activator, ref="noActivator")
# compute fit
fit_505 <- lm(RFU ~ time*activator, subset(data_505, time >= 10))
# predict trajectory
data_predict_505 <- data.frame(time=rep(seq(time_start, time_stop, 1), times=2),
                               activator=rep(c("noActivator", "100fM"), 
                                             each=length(seq(time_start, time_stop, 1))))
data_predict_505$RFU <- predict(fit_505, data_predict_505)
# pull regression coefficients
coef_505 <- data.frame(summary(fit_505)$coefficients)
# generate plot text
max_505 <- max(data_505$RFU, na.rm=T)
width_505 <- max(data_505$RFU, na.rm=T) - min(data_505$RFU, na.rm=T)
rate_noActivator_505 <- coef_505$Estimate[rownames(coef_505)=="time"]
rate_100fM_505 <- coef_505$Estimate[rownames(coef_505)=="time"] + 
  coef_505$Estimate[rownames(coef_505)=="time:activator100fM"]
pvalue_505 <- coef_505$Pr...t..[rownames(coef_505)=="time:activator100fM"]
plot_text_505 <- data.frame(x=10,
                            y=max_505 - c(1, 2.5, 4)*width_505/20,
                            text=c(paste("RNP-only rate:", signif(rate_noActivator_505, 3)),
                                   paste("100 fM rate:", signif(rate_100fM_505, 3)),
                                   paste("p =", signif(pvalue_505, 3))))
# generate plot
guide_plot_505 <- ggplot(data_predict_505, aes(x=time, y=RFU, col=activator)) +
  geom_line(size=2) + theme_classic(base_size=10) + 
  geom_point(data=data_505, shape=16, alpha=0.75,
             aes(x=time, y=RFU, col=activator, group=well_384)) +
  scale_color_manual(values=c("noActivator"="grey35", "100fM"="red")) +
  ggtitle("target region: 7720-7739",
          subtitle=paste(plot_text_505$text, collapse="\n")) + 
  xlab("time (min)") + ylab("RFU") + guides(col="none")

# aggregate plots ---------------------------------------------------------

figure_2B <- guide_plot_527 + guide_plot_505
ggsave(filename=file.path(figure_dir, "figure_2B.pdf"),
       plot=figure_2B,
       device="pdf", width=6.5, height=3, units="in")
