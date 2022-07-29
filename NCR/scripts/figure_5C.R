rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")
data_dir <- file.path(project_dir, "data", "supplementary_data")

# load plate_map -----------------------------------------------------------

plate_map_fname <- file.path(data_dir, 
                             "Anti-Tag Delta Omicron Pairwise Guide Screening (Modified)_20220719_platemap.csv")
plate_map <- read.csv(plate_map_fname)
plate_map <- plate_map[-1, ]
plate_map <- plate_map[, -1]
plate_map <- plate_map[, -ncol(plate_map)]
plate_map <- data.frame(plate_row = rep(LETTERS[seq(nrow(plate_map))], times=ncol(plate_map)),
                        plate_col = rep(seq(ncol(plate_map)), each=nrow(plate_map)),
                        sample = unlist(plate_map))
plate_map$plate_well <- with(plate_map, paste0(plate_row, plate_col))

plate_map$activator <- sapply(plate_map$plate_row,
                              function(x) {
                                ifelse(x %in% LETTERS[1:8], "Delta", "Omicron")
                              })

plate_map$sample <- sub("0fM", "0 fM", plate_map$sample)
plate_map$guide_id <- sub(" .*$", "", plate_map$sample)
plate_map$conc <- sapply(plate_map$sample, 
                         function(x) paste(strsplit(x, split=" ")[[1]][c(2,3)], collapse=" "))

plate_map <- subset(plate_map, grepl("^16", plate_map$guide_id))

# load guide features -----------------------------------------------------

guide_features_fname <- file.path(data_dir, 
                                  "Anti-Tag Delta Omicron Pairwise Guide Screening (Modified)_20220719_guides.csv")
guide_features <- read.csv(guide_features_fname)


plate_map$targeting <- sapply(plate_map$guide_id,
                              function(x) {
                                tmp_id <- paste0("NCR_", x)
                                ifelse(tmp_id %in% guide_features$delta.NCR.id, "Delta", 
                                       ifelse(tmp_id %in% guide_features$omicron.NCR.id,
                                              "Omicron", NA))
                              })
plate_map$variant <- sapply(seq(nrow(plate_map)),
                            function(x) {
                              tmp_id <- paste0("NCR_", plate_map$guide_id[x])
                              ifelse(plate_map$targeting[x] == "Delta",
                                     guide_features$variant[match(tmp_id, guide_features$delta.NCR.id)],
                                     guide_features$variant[match(tmp_id, guide_features$omicron.NCR.id)])
                            })
plate_map$tag <- sapply(seq(nrow(plate_map)),
                        function(x) {
                          tmp_id <- paste0("NCR_", plate_map$guide_id[x])
                          ifelse(plate_map$targeting[x]=="Delta",
                                 guide_features$delta.tag[match(tmp_id, guide_features$delta.NCR.id)],
                                 guide_features$omicron.tag[match(tmp_id, guide_features$omicron.NCR.id)])
                        })
plate_map$lineage <- sapply(seq(nrow(plate_map)),
                            function(x) {
                              tmp_id <- paste0("NCR_", plate_map$guide_id[x])
                              ifelse(plate_map$targeting[x]=="Delta",
                                     guide_features$twist.in.stock[match(tmp_id, guide_features$delta.NCR.id)],
                                     guide_features$twist.in.stock[match(tmp_id, guide_features$omicron.NCR.id)])
                            })

# load plate reader data --------------------------------------------------

plate_data_fname <- file.path(data_dir, 
                              "Anti-Tag Delta Omicron Pairwise Guide Screening (Modified)_20220719.xlsx")
plate_data <- readxl::read_xlsx(plate_data_fname,
                                range="A84:NW115")
plate_data$time <- round(plate_data$`Time [s]`/60) # convert from seconds to minutes

# parse data by guide
plate_data <- lapply(seq(nrow(plate_map)),
                     function(x) {
                       data.frame(time=plate_data$time,
                                  RFU=plate_data[[plate_map$plate_well[x]]],
                                  plate_cell=plate_map$plate_well[x],
                                  guide_id=plate_map$guide_id[x],
                                  conc=plate_map$conc[x],
                                  targeting=plate_map$targeting[x],
                                  variant=plate_map$variant[x],
                                  tag=plate_map$tag[x],
                                  activator=plate_map$activator[x])
                     })
plate_data <- do.call(rbind, plate_data)

# compute activator-dependent activity ------------------------------------

guide_activity <- lapply(unique(plate_map$guide_id),
                         function(x) {
                           tmp_data <- subset(plate_data,
                                              guide_id == x & time >= 10)
                           tmp_coef <- lapply(c("Delta", "Omicron"),
                                              function(y) {
                                                tmp_fit <- lm(RFU ~ time*conc,
                                                              subset(tmp_data, activator==y))
                                                tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                                tmp_coef$term <- rownames(tmp_coef)
                                                tmp_coef$guide_id <- x
                                                tmp_coef$activator <- y
                                                return(tmp_coef)
                                              })
                           tmp_coef <- do.call(rbind, tmp_coef)
                           return(tmp_coef)
                         })
guide_activity <- do.call(rbind, guide_activity)
guide_activity <- subset(guide_activity, grepl("time:", guide_activity$term))
guide_activity$variant <- plate_map$variant[match(guide_activity$guide_id,
                                                  plate_map$guide_id)]

detectable_variant <- sapply(split(guide_activity, guide_activity$variant), 
                             function(x) all(x$Pr...t.. < 0.05))
detectable_variant <- names(detectable_variant)[detectable_variant]

# same locus; both delta and omicron have an alternate allele at this position
detectable_variant <- detectable_variant[detectable_variant != "21618_G"] 


# compute guide_rates -----------------------------------------------------

guide_rates <- lapply(unique(subset(plate_map, variant %in% detectable_variant)$guide_id),
                      function(x) {
                        tmp_data <- subset(plate_data,
                                           conc == "167 fM" & guide_id == x & time >= 10)
                        tmp_coef <- lapply(c("Delta", "Omicron"),
                                           function(y) {
                                             tmp_fit <- lm(RFU ~ time,
                                                           data=subset(tmp_data,
                                                                       activator==y))
                                             tmp_coef <- data.frame(summary(tmp_fit)$coefficients)
                                             tmp_coef$term <- rownames(tmp_coef)
                                             tmp_coef$guide_id <- x
                                             tmp_coef$activator <- y
                                             return(tmp_coef)
                                           })
                        tmp_coef <- do.call(rbind, tmp_coef)
                        return(tmp_coef)
                      })
guide_rates <- do.call(rbind, guide_rates)
guide_rates <- subset(guide_rates, term=="time")
guide_rates <- plyr::join(guide_rates, 
                          unique(plate_map[, c("guide_id", "targeting", 
                                               "variant", "tag", "lineage")]), 
                          by="guide_id")

variant_label <- c("15240_T"="15240\nDelta: C\nOmicron: U",
                   "21618_T"="21618\nDelta: G\nOmicron: C",
                   "23202_A"="23202\nDelta: C\nOmicron: A")
guide_rates$variant_label <- variant_label[guide_rates$variant]

guide_labels <- c("1659"="Delta\nAAAG",
                  "1669"="Omicron\nAAAA",
                  "1660"="Delta\nAAAC",
                  "1670"="Omicron\nAAAG",
                  "1661"="Delta\nAAAG",
                  "1671"="Omicron\nAAAU")
guide_rates$guide_label <- guide_labels[guide_rates$guide_id]

figure_5C <- ggplot(guide_rates,
                    aes(x=guide_label, y=Estimate, fill=activator,
                        ymin=Estimate+qnorm(0.025)*Std..Error,
                        ymax=Estimate+qnorm(0.975)*Std..Error)) + 
  geom_col(position="dodge", width=0.8) + 
  geom_errorbar(position=position_dodge(0.8), width=0.5) + 
  theme_classic(base_size=8) + 
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) + 
  xlab("tag sequence (5' to 3')") + ylab("RFU/min") + 
  facet_wrap(~variant_label, scales="free_x") + 
  theme(axis.title.x=element_text(margin=margin(t=15, r=0, b=0, l=0)),
        strip.background=element_rect(fill="#FFA626")) + 
  scale_fill_manual(values=c("Delta"="#1B9E77", "Omicron"="#7570B3"))

ggsave(filename=file.path(figure_dir, "figure_5C.pdf"),
       plot=figure_5C,
       device="pdf", width=6.5, height=2, units="in")