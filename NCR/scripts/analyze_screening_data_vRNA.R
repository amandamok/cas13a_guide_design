rm(list=ls())

library(here)
library(ggplot2)

project_dir <- file.path(here(), "NCR")
ref_dir <- file.path(here(), "ref_data")
data_dir <- file.path(project_dir, "data")
trace_viral_dir <- file.path(project_dir, "figures", "traces", "viral")

guide_features <- read.csv(file.path(data_dir, "guide_features.csv"))

# functions ---------------------------------------------------------------

analyze_guide <- function(screening_data, guide, signal="RFU",
                          time_start=10, mixed_model=F) {
  # screening_data: data.frame
  # guide: character
  # signal: character; column in screening_data to be used for linear regression
  # model: formula to be used for lm() regression
  # time_start: timepoint to start regression
  # mixed_model: logical; whether to run linear regression with random effects for 384-well
  data_subset <- subset(screening_data, guide_id==guide)
  data_subset$activator <- relevel(data_subset$activator, ref="noActivator")
  time_stop <- signif(max(data_subset$time))
  if(!mixed_model) {
    # compute fit
    model_fit <- lm(formula(paste(signal, "~ time*activator")), 
                    data=subset(data_subset, time >= time_start))
    # predict trajectory
    model_fit_coef <- coef(model_fit)
    data_predict <- expand.grid(time=seq(time_start, time_stop, 1),
                                activator=c("noActivator", "100fM"))
    data_predict[[signal]] <- predict(model_fit, newdata=data_predict)
    # pull coefficients
    guide_coef <- data.frame(summary(model_fit)$coefficients, guide_id=guide,
                             coef=rownames(summary(model_fit)$coefficients), row.names=NULL)
    colnames(guide_coef)[colnames(guide_coef)=="Pr...t.."] <- "p.value"
    colnames(guide_coef)[colnames(guide_coef)=="Std..Error"] <- "Std.Error"
  } else {
    # compute fit
    model_fit <- nlme::lme(formula(paste(signal, "~ time*activator")),
                           random = ~ 1 + time | well_384, 
                           data = subset(data_subset, time >= time_start))
    # predict trajectory
    model_fit_coef <- summary(model_fit)$tTable[,1]
    data_predict_noActivator <- data.frame(time=seq(time_start, time_stop, 1),
                                           activator="noActivator")
    data_predict_noActivator[[signal]] <- model_fit_coef["(Intercept)"] + 
      data_predict_noActivator$time * model_fit_coef["time"]
    data_predict_100fM <- data.frame(time=seq(time_start, time_stop, 1),
                                     activator="100fM")
    data_predict_100fM[[signal]] <- model_fit_coef["(Intercept)"] + 
      model_fit_coef["activator100fM"] +
      data_predict_100fM$time * (model_fit_coef["time"] + model_fit_coef["time:activator100fM"])
    data_predict <- rbind(data_predict_noActivator, data_predict_100fM)
    # pull fixed effects
    guide_coef <- data.frame(summary(model_fit)$tTable, guide_id=guide,
                             coef=rownames(summary(model_fit)$varFix), row.names=NULL)
    colnames(guide_coef)[colnames(guide_coef)=="Value"] <- "Estimate"
    guide_coef <- guide_coef[, -which(colnames(guide_coef)=="DF")]
  }
  # make plot
  plot_max <- max(data_subset[[signal]], na.rm=T)
  plot_width <- max(data_subset[[signal]], na.rm=T) - min(data_subset[[signal]], na.rm=T)
  rate_noActivator <- model_fit_coef["time"]
  rate_100fM <- model_fit_coef["time"] + model_fit_coef["time:activator100fM"]
  plot_text <- data.frame(x=10,
                          y=plot_max - c(1:3)*plot_width/20,
                          text=c(paste("0 fM rate:", signif(rate_noActivator, 3)),
                                 paste("100 fM rate:", signif(rate_100fM, 3)),
                                 paste("p =", 
                                       signif(subset(guide_coef, coef=="time:activator100fM")$p.value, 3))))
  guide_plot <- ggplot(data_predict, 
                       aes_string(x="time", y=signal, col="activator")) +
    geom_line(size=3) + 
    geom_point(data=data_subset, 
               aes_string(x="time", y=signal, col="activator"), alpha=0.5) +
    geom_text(data=plot_text, aes(x=x, y=y, label=text), col="black", hjust=0) +
    theme_bw() + # scale_color_manual(values=c("red", "black")) +
    scale_color_manual(values=c("grey35", "red")) +
    ggtitle(paste0("NCR_", guide)) + xlab("time (min)") + ylab("RFU")
  return(list(plot=guide_plot, coef=guide_coef))
}

# generate mapping between 96well and 384well plates ----------------------

mapping_384_to_96 <- data.frame(matrix(NA, nrow=16, ncol=24))
rownames(mapping_384_to_96) <- LETTERS[1:16]
for(x in 1:8) {
  for(y in 1:12) {
    tmp_rows <- c(2*x-1, 2*x-1, 2*x)
    tmp_cols <- c(2*y-1, 2*y, 2*y)
    tmp_label <- paste0(LETTERS[x], y)
    for(z in 1:3) {
      mapping_384_to_96[tmp_rows[z], tmp_cols[z]] <- tmp_label
    }
  }
}
mapping_384_to_96 <- data.frame(well_96 = unlist(mapping_384_to_96),
                                well_384 = paste0(rep(LETTERS[1:16], times=24),
                                                  rep(1:24, each=16)))
mapping_384_to_96[is.na(mapping_384_to_96)] <- "empty"

# load platemaps ----------------------------------------------------------

expt_names <- c("GS1_1", "GS1_2", "GS2_1", "GS2_2", "GS3")
plate_maps <- openxlsx::read.xlsx(file.path(data_dir, "AC_Guide_Screening.xlsx"),
                                  sheet="Plate Setups")
plate_maps <- subset(plate_maps, X2 %in% LETTERS[1:8])[, 3:14]
plate_maps <- split(plate_maps, rep(expt_names, each=8))
plate_maps <- lapply(seq_along(plate_maps),
                     function(x) {
                       if(x %in% 1:4) {
                         data.frame(sample = paste0(sub("\\.0", "", 
                                                        unlist(plate_maps[[x]])), 
                                                    rep(c("_noActivator", "_100fM"), 
                                                        each=48)),
                                    well_96 = paste0(rep(LETTERS[1:8], times=12),
                                                     rep(1:12, each=8)))
                       } else {
                         tmp_map <- data.frame(sample = paste0(sub("\\.0", "", 
                                                                   unlist(plate_maps[[x]][1:4,])),
                                                               c(rep(rep(c("_noActivator", "_100fM"), 
                                                                         each=2),
                                                                     times=12))),
                                               well_96 = paste0(rep(LETTERS[1:4], times=12),
                                                                rep(1:12, each=4)))
                         return(subset(tmp_map, !grepl("NA", tmp_map$sample)))
                       }
                     })
names(plate_maps) <- expt_names

# load plate data ---------------------------------------------------------

screening_data_dir <- file.path(data_dir, "100fM_data")
data_fnames <- list.files(screening_data_dir)
plate_data <- lapply(data_fnames,
                     function(x) {
                       tmp_plate <- sub("_100fM_Plate", "", 
                                        sub("-", "_", 
                                            sub(".xlsx", "", x)))
                       row_indices <- 57:117
                       tmp_data <- openxlsx::read.xlsx(file.path(screening_data_dir, x), 
                                                       sheet=1, rows=row_indices)
                       tmp_data <- lapply(seq(nrow(tmp_data)),
                                          function(x) {
                                            data.frame(tmp_data[x, 1:3],
                                                       unlist(tmp_data[x, 4:387]),
                                                       colnames(tmp_data)[4:387],
                                                       row.names = NULL)
                                          })
                       tmp_data <- do.call(rbind, tmp_data)
                       colnames(tmp_data) <- c("cycle_number", "time", "temp", "RFU", "well_384")
                       tmp_data$well_96 <- mapping_384_to_96$well_96[match(tmp_data$well_384,
                                                                           mapping_384_to_96$well_384)]
                       tmp_platemap <- plate_maps[[tmp_plate]]
                       tmp_data$sample <- tmp_platemap$sample[match(tmp_data$well_96,
                                                                    tmp_platemap$well_96)]
                       tmp_data$plate <- tmp_plate
                       return(tmp_data)
                     })
plate_data <- do.call(rbind, plate_data)
plate_data$guide_id <- sapply(plate_data$sample,
                              function(x) {
                                ifelse(grepl("_", x),
                                       sub(" ", "_", sub("_.*", "", x)),
                                       "empty")
                              })

plate_data$activator <- factor(sapply(plate_data$sample,
                                      function(x) {
                                        ifelse(grepl("_", x),
                                               sub(".*_", "", x),
                                               "empty")
                                      }), 
                               levels=c("100fM", "noActivator", "empty"))
plate_data$cycle_number <- as.numeric(plate_data$cycle_number)
plate_data$RFU <- as.numeric(plate_data$RFU)
plate_data$time <- round(plate_data$time / 60) # convert from seconds to minutes
plate_data <- subset(plate_data, guide_id != "empty")

# remove outlier: C15 of 384 well plate (NCR_1324,  100fM) ----------------

plate_data <- subset(plate_data, well_384 != "C15")

# analyze data ------------------------------------------------------------

guide_ids <- unique(plate_data$guide_id)
guide_ids <- guide_ids[!(guide_ids %in% c("612_Control", "No_Protein"))]
all_results <- lapply(guide_ids,
                      function(x) {
                        analyze_guide(plate_data, guide=x, time_start=10, 
                                      signal="RFU", mixed_model=T)
                      })
names(all_results) <- guide_ids

# aggregate coefficients
all_coef <- lapply(all_results, function(x) x$coef)
all_coef <- do.call(rbind, all_coef)
all_coef <- cbind(all_coef,
                  guide_features[match(paste0("NCR_", all_coef$guide_id),
                                       guide_features$NCR.id),])
all_coef$plate <- plate_data$plate[match(all_coef$guide_id, plate_data$guide_id)]

guide_rate <- subset(all_coef, grepl(":", all_coef$coef))
noActivator_rate <- subset(all_coef, all_coef$coef=="time")

# save all traces ---------------------------------------------------------

all_plots <- lapply(seq_along(all_results),
                    function(x) {
                      guide_id <- paste0("NCR_", names(all_results)[x])
                      which_row <- match(guide_id, guide_features$NCR.id)
                      has_hairpin <- guide_features$has_crRNA_hairpin[which_row]
                      spacer_bp <- guide_features$crRNA_spacer_basepairs[which_row]
                      plot_subtitle <- paste(ifelse(has_hairpin, "has hairpin", "no hairpin"), "|",
                                             spacer_bp, "basepaired positions in spacer")
                      all_results[[x]]$plot + ggtitle(guide_id, subtitle=plot_subtitle)
                    })
names(all_plots) <- guide_ids

if(!dir.exists(trace_viral_dir)) { dir.create(trace_viral_dir) }
for(x in seq_along(all_plots)) {
  ggsave(filename=file.path(trace_viral_dir, 
                            paste0("NCR_", names(all_plots)[x], ".pdf")),
         plot=all_plots[[x]], device="pdf")
}

# write results to file ---------------------------------------------------

write.csv(guide_rate, file=file.path(data_dir, "guide_rate.csv"), row.names=F)
write.csv(guide_rate, file=file.path(data_dir, "noActivator_rate.csv"), row.names=F)
save(plate_data, file=file.path(data_dir, "vRNA_plate_data.Rda"))
