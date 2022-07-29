rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(latex2exp)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data", "supplementary_data")

# load data ---------------------------------------------------------------

data_fname <- file.path(data_dir, "EM_Cas13Lod_Rep3_x8_06152022_18966Zh_80of32gain.xlsx")

platemap <- data.frame(plate_row=rep(LETTERS[5:12], each=20),
                       plate_col=rep(c(3:22), times=8),
                       conc=rep(rep(c("no_protein", "no_activator",
                                      "125", "63", "31", "16", "8", "4", "2", "0"),
                                    each=2), times=8))
platemap$multiplex <- ifelse(platemap$plate_row %in% LETTERS[5:8], "8g", "32g")
platemap$well <- with(platemap, paste(plate_row, plate_col, sep=""))

raw_data <- readxl::read_xlsx(data_fname, range="A76:NW135")
raw_data <- reshape2::melt(raw_data, id.vars=1:3, variable.name="well")
raw_data$time <- raw_data$`Time [s]`/60
raw_data$value <- as.numeric(raw_data$value)
raw_data <- subset(raw_data, well %in% platemap$well)
raw_data <- droplevels(raw_data)
raw_data <- plyr::join(raw_data, platemap, by="well")

# exploratory analyses ----------------------------------------------------

all_rates <- lapply(unique(raw_data$well),
                    function(x) {
                      tmp_data <- subset(raw_data, well==x)
                      tmp_fit <- lm(value ~ time, subset(tmp_data, time > 10))
                      tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                      tmp_coef$term <- rownames(tmp_coef)
                      tmp_coef$well <- x
                      return(tmp_coef)
                    })
all_rates <- do.call(rbind, all_rates)
all_rates <- subset(all_rates, term=="time")
all_rates <- plyr::join(all_rates, platemap, by="well")

rates_first_half <- lapply(unique(raw_data$well),
                           function(x) {
                             tmp_data <- subset(raw_data, well==x)
                             tmp_fit <- lm(value ~ time,
                                           subset(tmp_data, 
                                                  time > 10 & time < 60))
                             tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                             tmp_coef$term <- rownames(tmp_coef)
                             tmp_coef$well <- x
                             return(tmp_coef)
                           })
rates_first_half <- do.call(rbind, rates_first_half)
rates_first_half <- subset(rates_first_half, term=="time")
rates_first_half <- plyr::join(rates_first_half, platemap, by="well")

rates_second_half <- lapply(unique(raw_data$well),
                            function(x) {
                              tmp_data <- subset(raw_data, well==x)
                              tmp_fit <- lm(value ~ time,
                                            subset(tmp_data, time > 60))
                              tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                              tmp_coef$term <- rownames(tmp_coef)
                              tmp_coef$well <- x
                              return(tmp_coef)
                            })
rates_second_half <- do.call(rbind, rates_second_half)
rates_second_half <- subset(rates_second_half, term=="time")
rates_second_half <- plyr::join(rates_second_half, platemap, by="well")

ggplot(raw_data, aes(x=time, y=value, group=well,
                     col=factor(conc, levels=unique(platemap$conc)))) + 
  geom_line(alpha=0.5) + theme_classic() + 
  xlab("time (min)") + ylab("RFU") + facet_grid( ~ multiplex) + labs(col="")

ggplot(all_rates, 
       aes(x=well, y=Estimate, fill=plate_row,
           ymin=Estimate + qnorm(0.025)*Std..Error,
           ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col() + geom_errorbar() + theme_classic() + 
  xlab("") + ylab("RFU/min") + guides(fill="none") + 
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

(ggplot(all_rates, 
        aes(x=well, y=Estimate, fill=plate_row,
            ymin=Estimate + qnorm(0.025)*Std..Error,
            ymax=Estimate + qnorm(0.975)*Std..Error)) + 
    geom_col() + geom_errorbar() + theme_classic() + 
    xlab("") + ylab("RFU/min") + guides(fill="none") + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
    ggtitle("time: 10-120 minutes")) / 
  (ggplot(rates_first_half, 
          aes(x=well, y=Estimate, fill=plate_row,
              ymin=Estimate + qnorm(0.025)*Std..Error,
              ymax=Estimate + qnorm(0.975)*Std..Error)) + 
     geom_col() + geom_errorbar() + theme_classic() + 
     xlab("") + ylab("RFU/min") + guides(fill="none") + 
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
     ggtitle("time: 10-60 minutes")) /
  ((ggplot(rates_second_half, 
           aes(x=well, y=Estimate, fill=plate_row,
               ymin=Estimate + qnorm(0.025)*Std..Error,
               ymax=Estimate + qnorm(0.975)*Std..Error)) + 
      geom_col() + geom_errorbar() + theme_classic() + 
      xlab("") + ylab("RFU/min") + guides(fill="none") + 
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
      ggtitle("time: 60-120 minutes")))

(ggplot(subset(raw_data,
               well %in% subset(all_rates, Estimate < 0)$well), 
        aes(x=time, y=value, col=well)) + 
    geom_line() + theme_classic() + guides(col="none") + 
    xlab("time (min)") + ylab("RFU") + 
    ggtitle("wells with negative slope")) /
  (ggplot(subset(all_rates, Estimate < 0),
          aes(x=well, y=Estimate, fill=well,
              ymin=Estimate + qnorm(0.025)*Std..Error,
              ymax=Estimate + qnorm(0.975)*Std..Error)) + 
     geom_col() + geom_errorbar() + theme_classic() + 
     xlab("") + ylab("RFU/min") + guides(fill="none") + 
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))

ggplot(subset(raw_data,
              well %in% subset(all_rates, Pr...t.. > 0.05)$well),
       aes(x=time, y=value, col=conc, group=well)) + 
  geom_line(alpha=0.5) + theme_classic() + 
  xlab("time (min)") + ylab("RFU")

# LOD rates: 10-120 minutes -----------------------------------------------

lod_rates <- lapply(unique(platemap$multiplex),
                    function(tmp_multiplex) {
                      tmp_rates <- lapply(unique(platemap$conc),
                                          function(tmp_conc) {
                                            tmp_data <- subset(raw_data, 
                                                               multiplex==tmp_multiplex & 
                                                                 conc==tmp_conc)
                                            tmp_fit <- lm(value ~ time, 
                                                          subset(tmp_data, time > 10))
                                            tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                            tmp_coef$term <- rownames(tmp_coef)
                                            tmp_coef$conc <- tmp_conc
                                            return(tmp_coef)
                                          })
                      tmp_rates <- do.call(rbind, tmp_rates)
                      tmp_rates$multiplex <- tmp_multiplex
                      return(tmp_rates)
                    })
lod_rates <- do.call(rbind, lod_rates)
lod_rates <- subset(lod_rates, term=="time")

lod_rates$conc <- factor(lod_rates$conc, levels=unique(platemap$conc))
lod_rates$multiplex <- factor(lod_rates$multiplex, levels=unique(platemap$multiplex))
lod_rates$detected <- sapply(seq(nrow(lod_rates)),
                             function(x) {
                               tmp_multiplex <- lod_rates$multiplex[x]
                               tmp_conc <- as.character(lod_rates$conc[x])
                               if(grepl("_", tmp_conc)) {
                                 return("")
                               } else {
                                 tmp_data <- rbind(subset(raw_data,
                                                          multiplex==tmp_multiplex &
                                                            conc==tmp_conc),
                                                   subset(raw_data,
                                                          multiplex==tmp_multiplex & 
                                                            conc=="no_activator"))
                                 tmp_data$conc <- factor(tmp_data$conc,
                                                         levels=c("no_activator", tmp_conc))
                                 tmp_fit <- lm(value ~ time*conc, 
                                               subset(tmp_data, time > 10))
                                 tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                 tmp_p <- tmp_coef$Pr...t..[grepl(":", rownames(tmp_coef))]
                                 tmp_rate <- tmp_coef$Estimate[grepl(":", rownames(tmp_coef))]
                                 return(ifelse(tmp_p/2 < 0.05 & tmp_rate > 0, "*", ""))
                               }
                             })

ggplot(lod_rates, 
       aes(x=conc, y=Estimate, fill=multiplex,
           ymin=Estimate + qnorm(0.025)*Std..Error,
           ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col(position=position_dodge()) + 
  geom_errorbar(position=position_dodge(width=0.9), width=0.5) + 
  theme_classic() + scale_fill_manual(values=c("8g"="grey35", "32g"="red")) +  
  xlab("") + ylab("RFU/min") + labs(fill="") + 
  geom_text(aes(y=Estimate + qnorm(0.975)*Std..Error + 5,
                label=detected), position=position_dodge(width=0.9))

# LOD rates: 10-60 minutes ------------------------------------------------

lod_rates_firsthalf <- lapply(unique(platemap$multiplex),
                              function(tmp_multiplex) {
                                tmp_rates <- lapply(unique(platemap$conc),
                                                    function(tmp_conc) {
                                                      tmp_data <- subset(raw_data, 
                                                                         multiplex==tmp_multiplex & 
                                                                           conc==tmp_conc)
                                                      tmp_fit <- lm(value ~ time, 
                                                                    subset(tmp_data, 
                                                                           time > 10 & time < 60))
                                                      tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                                      tmp_coef$term <- rownames(tmp_coef)
                                                      tmp_coef$conc <- tmp_conc
                                                      return(tmp_coef)
                                                    })
                                tmp_rates <- do.call(rbind, tmp_rates)
                                tmp_rates$multiplex <- tmp_multiplex
                                return(tmp_rates)
                              })
lod_rates_firsthalf <- do.call(rbind, lod_rates_firsthalf)
lod_rates_firsthalf <- subset(lod_rates_firsthalf, term=="time")

lod_rates_firsthalf$conc <- factor(lod_rates_firsthalf$conc, levels=unique(platemap$conc))
lod_rates_firsthalf$multiplex <- factor(lod_rates_firsthalf$multiplex, levels=unique(platemap$multiplex))
lod_rates_firsthalf$detected <- sapply(seq(nrow(lod_rates_firsthalf)),
                                       function(x) {
                                         tmp_multiplex <- lod_rates_firsthalf$multiplex[x]
                                         tmp_conc <- as.character(lod_rates_firsthalf$conc[x])
                                         if(grepl("_", tmp_conc)) {
                                           return("")
                                         } else {
                                           tmp_data <- rbind(subset(raw_data,
                                                                    multiplex==tmp_multiplex &
                                                                      conc==tmp_conc),
                                                             subset(raw_data,
                                                                    multiplex==tmp_multiplex & 
                                                                      conc=="no_activator"))
                                           tmp_data$conc <- factor(tmp_data$conc,
                                                                   levels=c("no_activator", tmp_conc))
                                           tmp_fit <- lm(value ~ time*conc, 
                                                         subset(tmp_data, 
                                                                time > 10 & time < 60))
                                           tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                           tmp_p <- tmp_coef$Pr...t..[grepl(":", rownames(tmp_coef))]
                                           tmp_rate <- tmp_coef$Estimate[grepl(":", rownames(tmp_coef))]
                                           return(ifelse(tmp_p/2 < 0.05 & tmp_rate > 0, "*", ""))
                                         }
                                       })

ggplot(lod_rates_firsthalf, 
       aes(x=conc, y=Estimate, fill=multiplex,
           ymin=Estimate + qnorm(0.025)*Std..Error,
           ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col(position=position_dodge()) + 
  geom_errorbar(position=position_dodge(width=0.9), width=0.5) + 
  theme_classic() + scale_fill_manual(values=c("8g"="grey35", "32g"="red")) +  
  xlab("") + ylab("RFU/min") + labs(fill="") + 
  geom_text(aes(y=Estimate + qnorm(0.975)*Std..Error + 5,
                label=detected), position=position_dodge(width=0.9))

# LOD rates: 60-120 minutes -----------------------------------------------

lod_rates_secondhalf <- lapply(unique(platemap$multiplex),
                               function(tmp_multiplex) {
                                 tmp_rates <- lapply(unique(platemap$conc),
                                                     function(tmp_conc) {
                                                       tmp_data <- subset(raw_data, 
                                                                          multiplex==tmp_multiplex & 
                                                                            conc==tmp_conc)
                                                       tmp_fit <- lm(value ~ time, 
                                                                     subset(tmp_data, time > 60))
                                                       tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                                       tmp_coef$term <- rownames(tmp_coef)
                                                       tmp_coef$conc <- tmp_conc
                                                       return(tmp_coef)
                                                     })
                                 tmp_rates <- do.call(rbind, tmp_rates)
                                 tmp_rates$multiplex <- tmp_multiplex
                                 return(tmp_rates)
                               })
lod_rates_secondhalf <- do.call(rbind, lod_rates_secondhalf)
lod_rates_secondhalf <- subset(lod_rates_secondhalf, term=="time")

lod_rates_secondhalf$conc <- factor(lod_rates_secondhalf$conc, levels=unique(platemap$conc))
lod_rates_secondhalf$multiplex <- factor(lod_rates_secondhalf$multiplex, levels=unique(platemap$multiplex))
lod_rates_secondhalf$detected <- sapply(seq(nrow(lod_rates_secondhalf)),
                                        function(x) {
                                          tmp_multiplex <- lod_rates_secondhalf$multiplex[x]
                                          tmp_conc <- as.character(lod_rates_secondhalf$conc[x])
                                          if(grepl("_", tmp_conc)) {
                                            return("")
                                          } else {
                                            tmp_data <- rbind(subset(raw_data,
                                                                     multiplex==tmp_multiplex &
                                                                       conc==tmp_conc),
                                                              subset(raw_data,
                                                                     multiplex==tmp_multiplex & 
                                                                       conc=="no_activator"))
                                            tmp_data$conc <- factor(tmp_data$conc,
                                                                    levels=c("no_activator", tmp_conc))
                                            tmp_fit <- lm(value ~ time*conc, 
                                                          subset(tmp_data, time > 60))
                                            tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                            tmp_p <- tmp_coef$Pr...t..[grepl(":", rownames(tmp_coef))]
                                            tmp_rate <- tmp_coef$Estimate[grepl(":", rownames(tmp_coef))]
                                            return(ifelse(tmp_p/2 < 0.05 & tmp_rate > 0, "*", ""))
                                          }
                                        })

ggplot(lod_rates_secondhalf, 
       aes(x=conc, y=Estimate, fill=multiplex,
           ymin=Estimate + qnorm(0.025)*Std..Error,
           ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col(position=position_dodge()) + 
  geom_errorbar(position=position_dodge(width=0.9), width=0.5) + 
  theme_classic() + scale_fill_manual(values=c("8g"="grey35", "32g"="red")) +  
  xlab("") + ylab("RFU/min") + labs(fill="") + 
  geom_text(aes(y=Estimate + qnorm(0.975)*Std..Error + 5,
                label=detected), position=position_dodge(width=0.9))



(ggplot(lod_rates, 
        aes(x=conc, y=Estimate, fill=multiplex,
            ymin=Estimate + qnorm(0.025)*Std..Error,
            ymax=Estimate + qnorm(0.975)*Std..Error)) + 
    geom_col(position=position_dodge()) + 
    geom_errorbar(position=position_dodge(width=0.9), width=0.5) + 
    theme_classic() + scale_fill_manual(values=c("8g"="grey35", "32g"="red")) +  
    xlab("") + ylab("RFU/min") + labs(fill="") + 
    geom_text(aes(y=Estimate + qnorm(0.975)*Std..Error + 5,
                  label=detected), position=position_dodge(width=0.9)) + 
    ggtitle("time: 10-120 minutes")) /
  (ggplot(lod_rates_firsthalf, 
          aes(x=conc, y=Estimate, fill=multiplex,
              ymin=Estimate + qnorm(0.025)*Std..Error,
              ymax=Estimate + qnorm(0.975)*Std..Error)) + 
     geom_col(position=position_dodge()) + 
     geom_errorbar(position=position_dodge(width=0.9), width=0.5) + 
     theme_classic() + scale_fill_manual(values=c("8g"="grey35", "32g"="red")) +  
     xlab("") + ylab("RFU/min") + labs(fill="") + 
     geom_text(aes(y=Estimate + qnorm(0.975)*Std..Error + 5,
                   label=detected), position=position_dodge(width=0.9)) + 
     ggtitle("time: 10-60 minutes")) / 
  (ggplot(lod_rates_secondhalf, 
          aes(x=conc, y=Estimate, fill=multiplex,
              ymin=Estimate + qnorm(0.025)*Std..Error,
              ymax=Estimate + qnorm(0.975)*Std..Error)) + 
     geom_col(position=position_dodge()) + 
     geom_errorbar(position=position_dodge(width=0.9), width=0.5) + 
     theme_classic() + scale_fill_manual(values=c("8g"="grey35", "32g"="red")) +  
     xlab("") + ylab("RFU/min") + labs(fill="") + 
     geom_text(aes(y=Estimate + qnorm(0.975)*Std..Error + 5,
                   label=detected), position=position_dodge(width=0.9)) + 
     ggtitle("time: 60-120 minutes"))
