rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data")

time_start <- 10
time_stop <- 120

# load plate data ---------------------------------------------------------

load(file.path(data_dir, "vRNA_plate_data.Rda"))

# good guide: NCR_527 -----------------------------------------------------

data_527 <- subset(plate_data, guide_id=="527")
data_527$activator <- droplevels(relevel(data_527$activator, ref="noActivator"))
# compute fit
fit_527 <- nlme::lme(formula(RFU ~ time*activator),
                     random = ~ 1 + time | well_384, 
                     data = subset(data_527, time >= time_start))
# pull regression coefficients
coef_527 <- data.frame(summary(fit_527)$tTable)
# predict trajectory
data_predict_527 <- data.frame(time=rep(seq(time_start, time_stop, 1), times=2),
                               activator=rep(c("noActivator", "100fM"), 
                                             each=length(seq(time_start, time_stop, 1))))
data_predict_527$RFU <- coef_527$Value[rownames(coef_527)=="(Intercept)"] + 
  data_predict_527$time * coef_527$Value[rownames(coef_527)=="time"]
data_predict_527$RFU[data_predict_527$activator=="100fM"] <- 
  data_predict_527$RFU[data_predict_527$activator=="100fM"] + 
  data_predict_527$time[data_predict_527$activator=="100fM"] *
  coef_527$Value[rownames(coef_527)=="time:activator100fM"]

# generate plot
guide_plot_527 <- ggplot(data_predict_527, aes(x=time, y=RFU, col=activator)) +
  geom_line(size=2, alpha=0.5) + theme_classic(base_size=8) + 
  geom_line(data=data_527, aes(x=time, y=RFU, col=activator, group=well_384)) +
  scale_color_manual(values=c("noActivator"="grey35", "100fM"="red")) +
  ggtitle("crRNA 527") +
  xlab("time (min)") + ylab("RFU") + guides(col="none")

barplot_527 <- ggplot(data.frame(rate=c(coef_527$Value[rownames(coef_527)=="time"],
                                        sum(coef_527$Value[rownames(coef_527) %in% 
                                                                c("time", "time:activator100fM")])),
                                 std_error=c(coef_527$Std.Error[rownames(coef_527)=="time"],
                                             coef_527$Std.Error[rownames(coef_527)=="time:activator100fM"]),
                                 sample=factor(c("RNP\nonly", "100\nfM"),
                                               levels=c("RNP\nonly", "100\nfM"))),
                      aes(x=sample, y=rate, fill=sample,
                          ymin=rate+qnorm(0.025)*std_error,
                          ymax=rate+qnorm(0.975)*std_error)) + 
  geom_col() + geom_errorbar(width=0.5) + theme_classic(base_size=8) + 
  scale_fill_manual(values=c("RNP\nonly"="grey35", "100\nfM"="red")) + 
  guides(fill="none") + xlab("") + ylab("RFU/min")

# bad guide: NCR_553 ------------------------------------------------------

data_553 <- subset(plate_data, guide_id=="553")
data_553$activator <- droplevels(relevel(data_553$activator, ref="noActivator"))
# compute fit
fit_553 <- nlme::lme(formula(RFU ~ time*activator),
                     random = ~ 1 + time | well_384, 
                     data = subset(data_553, time >= time_start))
# pull regression coefficients
coef_553 <- data.frame(summary(fit_553)$tTable)
# predict trajectory
data_predict_553 <- data.frame(time=rep(seq(time_start, time_stop, 1), times=2),
                               activator=rep(c("noActivator", "100fM"), 
                                             each=length(seq(time_start, time_stop, 1))))
data_predict_553$RFU <- coef_553$Value[rownames(coef_553)=="(Intercept)"] + 
  data_predict_553$time * coef_553$Value[rownames(coef_553)=="time"]
data_predict_553$RFU[data_predict_553$activator=="100fM"] <- 
  data_predict_553$RFU[data_predict_553$activator=="100fM"] + 
  data_predict_553$time[data_predict_553$activator=="100fM"] *
  coef_553$Value[rownames(coef_553)=="time:activator100fM"]

# generate plot
guide_plot_553 <- ggplot(data_predict_553, aes(x=time, y=RFU, col=activator)) +
  geom_line(size=2, alpha=0.5) + theme_classic(base_size=8) + 
  geom_line(data=data_553, aes(x=time, y=RFU, col=activator, group=well_384)) +
  scale_color_manual(values=c("noActivator"="grey35", "100fM"="red")) +
  ggtitle("crRNA 553") +
  xlab("time (min)") + ylab("RFU") + guides(col="none")

barplot_553 <- ggplot(data.frame(rate=c(coef_553$Value[rownames(coef_553)=="time"],
                                        sum(coef_553$Value[rownames(coef_553) %in% 
                                                             c("time", "time:activator100fM")])),
                                 std_error=c(coef_553$Std.Error[rownames(coef_553)=="time"],
                                             coef_553$Std.Error[rownames(coef_553)=="time:activator100fM"]),
                                 sample=factor(c("RNP\nonly", "100\nfM"),
                                               levels=c("RNP\nonly", "100\nfM"))),
                      aes(x=sample, y=rate, fill=sample,
                          ymin=rate+qnorm(0.025)*std_error,
                          ymax=rate+qnorm(0.975)*std_error)) + 
  geom_col() + geom_errorbar(width=0.5) + theme_classic(base_size=8) + 
  scale_fill_manual(values=c("RNP\nonly"="grey35", "100\nfM"="red")) + 
  guides(fill="none") + xlab("") + ylab("RFU/min")

# aggregate plots ---------------------------------------------------------

figure_2B <- guide_plot_527 + barplot_527 + guide_plot_553 + barplot_553 + 
  plot_layout(nrow=1, widths=c(3, 1, 3, 1))
ggsave(filename=file.path(figure_dir, "figure_2B.pdf"),
       plot=figure_2B,
       device="pdf", width=6.5, height=2, units="in")
