rm(list=ls())

library(here)
library(ggplot2)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
trace_gblock_dir <- file.path(project_dir, "figures", "traces", "gblock")

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
plate_maps <- xlsx::read.xlsx(file.path(data_dir, "AC_Guide_Screening.xlsx"),
                              sheetName="Plate Setups")
plate_maps <- subset(plate_maps, NA. %in% LETTERS[1:8])[, 3:14]
plate_maps <- split(plate_maps, rep(expt_names, each=8))
plate_maps <- lapply(seq_along(plate_maps),
                     function(x) {
                       if(x %in% 1:4) {
                         data.frame(sample = paste0(unlist(plate_maps[[x]]), 
                                                    rep(c("_noActivator", "_100fM"), each=48)),
                                    well_96 = paste0(rep(LETTERS[1:8], times=12),
                                                     rep(1:12, each=8)))
                       } else {
                         tmp_map <- data.frame(sample = paste0(unlist(plate_maps[[x]][1:4,]),
                                                               c(rep(rep(c("_noActivator", "_100fM"), 
                                                                         each=2),
                                                                     times=12))),
                                               well_96 = paste0(rep(LETTERS[1:4], times=12),
                                                                rep(1:12, each=4)))
                         return(subset(tmp_map, !grepl("NA", tmp_map$sample)))
                       }
                     })
names(plate_maps) <- expt_names

# load plate 1-2 platemap
gblock_platemap_1_2 <- xlsx::read.xlsx(file.path(data_dir, 
                                                 "Plate1_2_gBlock_384w_replicate_platemap.xlsx"),
                                       sheetIndex=1)
gblock_platemap_1_2 <- subset(gblock_platemap_1_2,
                              subset=value_type %in% LETTERS[1:16])[, -1]
gblock_platemap_1_2_activator <- data.frame(matrix("_noActivator", 
                                                   nrow=nrow(gblock_platemap_1_2),
                                                   ncol=ncol(gblock_platemap_1_2)))
gblock_platemap_1_2_activator[, 13:24] <- "_100fM"
gblock_platemap_1_2_activator[c(2, 4, 6), 9] <- "_100fM"
gblock_platemap_1_2_activator[c(2, 4, 6), 11] <- "_100fM"
gblock_platemap_1_2 <- data.frame(sample = paste0(unlist(gblock_platemap_1_2),
                                                  unlist(gblock_platemap_1_2_activator)),
                                  well_384 = paste0(rep(LETTERS[1:16], times=24),
                                                    rep(1:24, each=16)))
gblock_platemap_1_2$sample <- sub(" ", "_", gblock_platemap_1_2$sample)

# load plate 2-2 platemap
gblock_platemap_2_2 <- xlsx::read.xlsx(file.path(data_dir, 
                                                 "Plate2_2_gBlock_384w_replicate_platemap.xlsx"),
                                       sheetIndex=1)
gblock_platemap_2_2 <- subset(gblock_platemap_2_2,
                              subset=value_type %in% LETTERS[1:16])[, -1]
gblock_platemap_2_2_activator <- data.frame(matrix("_noActivator", 
                                                   nrow=nrow(gblock_platemap_2_2),
                                                   ncol=ncol(gblock_platemap_2_2)))
gblock_platemap_2_2_activator[, 13:24] <- "_100fM"
gblock_platemap_2_2_activator[c(2, 4, 6), 9] <- "_100fM"
gblock_platemap_2_2_activator[c(2, 4, 6), 11] <- "_100fM"
gblock_platemap_2_2 <- data.frame(sample = paste0(unlist(gblock_platemap_2_2),
                                                  unlist(gblock_platemap_2_2_activator)),
                                  well_384 = paste0(rep(LETTERS[1:16], times=24),
                                                    rep(1:24, each=16)))
gblock_platemap_2_2$sample <- sub(" ", "_", gblock_platemap_2_2$sample)

# load data ---------------------------------------------------------------

gblock_data_dir <- file.path(data_dir, "100fM_gblock_data")
gblock_data_fnames <- list.files(gblock_data_dir)
gblock_data <- lapply(gblock_data_fnames,
                      function(x) {
                        tmp_plate <- sub("_100fM_Plate", "",
                                         sub("-", "_",
                                             sub(".xlsx", "", x)))
                        row_indices <- 55:115
                        tmp_data <- xlsx::read.xlsx(file.path(gblock_data_dir, x),
                                                    sheetIndex=1, rowIndex=row_indices, header=T)
                        tmp_data <- lapply(seq(nrow(tmp_data)),
                                           function(x) {
                                             data.frame(tmp_data[x, 1:3],
                                                        unlist(tmp_data[x, 4:387]),
                                                        colnames(tmp_data)[4:387],
                                                        row.names=NULL)
                                           })
                        tmp_data <- do.call(rbind, tmp_data)
                        colnames(tmp_data) <- c("cycle_number", "time", "temp", "RFU", "well_384")
                        if(tmp_plate %in% paste0("gBlock", c("1_1", "2_1"))) {
                          well_96 <- mapping_384_to_96$well_96[match(tmp_data$well_384,
                                                                     mapping_384_to_96$well_384)]
                          tmp_platemap <- plate_maps[[sub("gBlock", "GS", tmp_plate)]]
                          tmp_data$sample <- tmp_platemap$sample[match(well_96,
                                                                       tmp_platemap$well_96)]
                        } else {
                          if(tmp_plate=="gBlock1_2") {
                            tmp_platemap <- gblock_platemap_1_2
                          } else {
                            tmp_platemap <- gblock_platemap_2_2
                          }
                          tmp_data$sample <- tmp_platemap$sample[match(tmp_data$well_384, tmp_platemap$well_384)]
                        }
                        tmp_data$plate <- tmp_plate
                        return(tmp_data)
                      })
gblock_data <- do.call(rbind, gblock_data)
gblock_data$sample[grepl("NA", gblock_data$sample)] <- "empty"
gblock_data$guide_id <- sapply(gblock_data$sample,
                               function(x) {
                                 ifelse(grepl("_", x),
                                        sub(" ", "_", sub("_.*", "", x)),
                                        "empty")
                               })
gblock_data$guide_id[gblock_data$guide_id=="612"] <- "612_Control"
gblock_data$guide_id[gblock_data$guide_id=="No"] <- "No_Protein"
gblock_data$activator <- factor(sapply(gblock_data$sample,
                                       function(x) {
                                         ifelse(grepl("_", x),
                                                sub(".*_", "", x),
                                                "empty")
                                       }), 
                                levels=c("100fM", "noActivator", "empty"))
gblock_data$cycle_number <- as.numeric(gblock_data$cycle_number)
gblock_data$RFU <- as.numeric(gblock_data$RFU)
gblock_data$time <- round(gblock_data$time/60) # convert from seconds to minutes
gblock_data <- subset(gblock_data, !(gblock_data$sample == "empty"))
gblock_data$activator <- droplevels(gblock_data$activator)
gblock_data <- subset(gblock_data, !is.na(RFU))

# analyze data ------------------------------------------------------------

gblock_guide_ids <- unique(gblock_data$guide_id)
gblock_guide_ids <- gblock_guide_ids[!(gblock_guide_ids %in% c("612_Control", "No_Protein"))]
gblock_results <- lapply(gblock_guide_ids,
                         function(x) {
                           analyze_guide(gblock_data, guide=x, time_start=10, 
                                         signal="RFU", mixed_model=F)
                         })
names(gblock_results) <- gblock_guide_ids

# aggregate coefficients
gblock_coef <- lapply(gblock_results, function(x) x$coef)
gblock_coef <- do.call(rbind, gblock_coef)
gblock_rate <- subset(gblock_coef, grepl(":", gblock_coef$coef))

# save all traces ---------------------------------------------------------

gblock_plots <- lapply(seq_along(gblock_results),
                       function(x) {
                         guide_id <- paste0("NCR_", names(gblock_results)[x])
                         which_row <- match(guide_id, guide_features$NCR.id)
                         has_hairpin <- guide_features$has_crRNA_hairpin[which_row]
                         spacer_bp <- guide_features$crRNA_spacer_basepairs[which_row]
                         plot_subtitle <- paste(ifelse(has_hairpin, "has hairpin", "no hairpin"), "|",
                                                spacer_bp, "basepaired positions in spacer")
                         gblock_results[[x]]$plot + ggtitle(guide_id, subtitle=plot_subtitle)
                       })
names(gblock_plots) <- gblock_guide_ids

if(!dir.exists(trace_gblock_dir)) { dir.create(trace_gblock_dir) }
for(x in seq_along(gblock_plots)) {
  ggsave(filename=file.path(trace_gblock_dir, 
                            paste0("NCR_", names(gblock_plots)[x], ".pdf")),
         plot=gblock_plots[[x]], device="pdf")
}

# write results to file ---------------------------------------------------

write.csv(gblock_rate, file=file.path(data_dir, "gblock_rate.csv"), row.names=F)
save(gblock_data, file=file.path(data_dir, "gblock_plate_data.Rda"))
