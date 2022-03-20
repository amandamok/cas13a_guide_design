rm(list=ls())

library(here)
library(ggplot2)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data", "supplementary_data")
figure_dir <- file.path(project_dir, "figures")

# load platemap -----------------------------------------------------------

lod_conc <- c(250, 125, 63, 31, 16, 8, 4, 0)
primary_lod_platemap <- data.frame(plate_row=rep(LETTERS[c(5, 7, 10, 11)], 
                                                 each=2*length(lod_conc)),
                                   plate_col=rep(5:20, times=4),
                                   multiplex=c(rep(8, 2*length(lod_conc)),
                                               rep(c(NA, 8), times=length(lod_conc)),
                                               rep(32, 2*length(lod_conc)),
                                               rep(c(32, NA), times=length(lod_conc))),
                                   conc=rep(rep(lod_conc, each=2), times=4))
# all replicates
primary_lod_platemap <- data.frame(plate_row=rep(LETTERS[5:12],
                                                 each=2*length(lod_conc)),
                                   plate_col=rep(5:20, times=8),
                                   multiplex=c(rep(8, 4*2*length(lod_conc)),
                                               rep(32, 4*2*length(lod_conc))),
                                   conc=rep(rep(lod_conc, each=2), times=8))

primary_lod_platemap$plate_well <- with(primary_lod_platemap, 
                                        paste0(plate_row, plate_col))
primary_lod_platemap <- subset(primary_lod_platemap, !is.na(multiplex))

# load data ---------------------------------------------------------------

primary_lod_data <- xlsx::read.xlsx(file.path(data_dir, "primary_LOD_20220204.xlsx"),
                                    sheetIndex=1, startRow=76, endRow=136, header=T)
primary_lod_data <- lapply(seq(nrow(primary_lod_platemap)),
                           function(x) {
                             tmp_well <- primary_lod_platemap$plate_well[x]
                             data.frame(time=primary_lod_data$Time..s.,
                                        RFU=as.numeric(primary_lod_data[[tmp_well]]),
                                        well=tmp_well,
                                        multiplex=primary_lod_platemap$multiplex[x],
                                        conc=primary_lod_platemap$conc[x])
                           })
primary_lod_data <- do.call(rbind, primary_lod_data)
primary_lod_data$time <- round(primary_lod_data$time / 60)

# compute slopes ----------------------------------------------------------

primary_lod_slopes <- lapply(unique(primary_lod_data$multiplex),
                             function(x) {
                               tmp_slopes <- lapply(unique(primary_lod_data$conc),
                                                    function(y) {
                                                      tmp_data <- subset(primary_lod_data,
                                                                         multiplex==x & 
                                                                           conc==y & 
                                                                           time >= 10)
                                                      tmp_model <- lm(RFU ~ time, tmp_data)
                                                      tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                      tmp_coef$coef <- rownames(tmp_coef)
                                                      tmp_coef$multiplex <- x
                                                      tmp_coef$conc <- y
                                                      return(tmp_coef)
                                                    })
                               tmp_slopes <- do.call(rbind, tmp_slopes)
                               return(tmp_slopes)
                             })
primary_lod_slopes <- do.call(rbind, primary_lod_slopes)
primary_lod_slopes <- subset(primary_lod_slopes, coef=="time")
primary_lod_slopes$conc <- factor(primary_lod_slopes$conc, 
                                  levels=unique(primary_lod_slopes$conc))


# compute individual slopes -----------------------------------------------

primary_lod_individual_slopes <- lapply(primary_lod_platemap$plate_well,
                                        function(x) {
                                          tmp_data <- subset(primary_lod_data, well==x)
                                          tmp_model <- lm(RFU ~ time, tmp_data)
                                          tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                          tmp_coef$coef <- rownames(tmp_coef)
                                          tmp_coef$plate_well <- x
                                          return(tmp_coef)
                                        })
primary_lod_individual_slopes <- do.call(rbind, primary_lod_individual_slopes)
primary_lod_individual_slopes <- subset(primary_lod_individual_slopes, coef=="time")
primary_lod_individual_slopes <- dplyr::left_join(primary_lod_individual_slopes,
                                                  primary_lod_platemap,
                                                  by="plate_well")
ggplot(primary_lod_individual_slopes, 
       aes(x=factor(conc), y=Estimate,
           ymin=Estimate+qnorm(0.025)*Std..Error,
           ymax=Estimate+qnorm(0.975)*Std..Error)) + 
  geom_col(position="dodge2") + geom_errorbar(position="dodge2") + 
  facet_grid(multiplex~., scales="free_y") + theme_classic() + 
  xlab("copies/uL") + ylab("RFU/min") + 
  geom_hline(yintercept=0)

# compare slopes between 8-plex and 32-plex -------------------------------

primary_lod_compare_slopes <- lapply(c(8, 32),
                                     function(x) {
                                       tmp <- lapply(lod_conc[-length(lod_conc)],
                                                     function(y) {
                                                       tmp_data <- subset(primary_lod_data,
                                                                          multiplex==x & 
                                                                            conc %in% c(0, y) &
                                                                            time >= 10)
                                                       tmp_data$conc <- factor(tmp_data$conc, 
                                                                               levels=c(0, y))
                                                       tmp_model <- lm(RFU ~ time*conc, tmp_data)
                                                       tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                       tmp_coef$coef <- rownames(tmp_coef)
                                                       tmp_coef$multiplex <- x
                                                       tmp_coef$conc <- y
                                                       return(tmp_coef)
                                                     })
                                       return(do.call(rbind, tmp))
                                     })
primary_lod_compare_slopes <- do.call(rbind, primary_lod_compare_slopes)
primary_lod_compare_slopes <- subset(primary_lod_compare_slopes, 
                                     grepl("time:", primary_lod_compare_slopes$coef))

# add pvalue annotation
primary_lod_slopes$lod_pvalue <- primary_lod_compare_slopes$Pr...t..[prodlim::row.match(primary_lod_slopes[, c("multiplex", "conc")],
                                                                                        primary_lod_compare_slopes[, c("multiplex", "conc")])]
primary_lod_slopes$lod_pvalue_label <- cut(primary_lod_slopes$lod_pvalue, 
                                           breaks=c(0, 0.001, 0.01, 0.05, 1),
                                           labels=c("p<0.001", 
                                                    "p<0.01", 
                                                    "p<0.05", 
                                                    "p>=0.05"))

# generate plot -----------------------------------------------------------

ggplot(primary_lod_slopes, 
       aes(x=conc, y=Estimate, 
           ymin=Estimate+qnorm(0.025)*Std..Error,
           ymax=Estimate+qnorm(0.975)*Std..Error, 
           fill=lod_pvalue_label)) + 
  geom_col() + geom_errorbar() + theme_classic() + 
  facet_grid(~multiplex) + xlab("") + ylab("RFU/min") + labs(fill="")

