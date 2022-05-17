rm(list=ls())

library(here)

data_dir <- file.path(here(), "NCR", "data", "supplementary_data")
platemap_fname <-"Anti-Tag experiment Platemap 05.12.22.xlsx"
data_fname <- "Anti-Tag Cas13 primary 051222 Guide_Screening (Modified)_20220512.xlsx"
mappings_fname <- "Anti-Tag experiment mappings.csv"

complementary_nt <- setNames(c("A", "T", "G", "C"), c("U", "A", "C", "G"))

# load platemap -----------------------------------------------------------

# read in file
platemap <- readxl::read_xlsx(file.path(data_dir, platemap_fname),
                              range="B2:Y12")
platemap <- data.frame(label=unlist(platemap),
                       plate_row=rep(LETTERS[1:nrow(platemap)], times=ncol(platemap)),
                       plate_col=rep(seq(ncol(platemap)), each=nrow(platemap)),
                       row.names=NULL)
platemap <- subset(platemap, !is.na(platemap$label))

# add annotations
platemap$label <- sub(" fM", "fM", platemap$label)
platemap$label <- sub("No protein", "NoProtein", platemap$label)
platemap$plate_cell <- with(platemap, paste0(plate_row, plate_col))
platemap$guide_id <- paste0("NCR_", sub(" .*$", "", platemap$label))
platemap$guide_id[platemap$guide_id == "NCR_NoProtein"] <- "NoProtein"
platemap$conc <- sub("^.* ", "", platemap$label)
platemap$conc[!(platemap$conc %in% c("0fM", "167fM"))] <- NA

# load mappings -----------------------------------------------------------

# read in file
mappings <- read.csv(file.path(data_dir, mappings_fname))

# add annotations
platemap$group <- mappings$NCR.id[match(platemap$guide_id,
                                        mappings$new.NCR.id)]
platemap$group[is.na(platemap$group)] <- platemap$guide_id[is.na(platemap$group)]
platemap$tag <- mappings$new_tag_pos1[match(platemap$guide_id, 
                                            mappings$new.NCR.id)]
platemap$tag[is.na(platemap$tag) & !(platemap$guide_id == "NoProtein")] <- "C"
platemap$antitag <- sapply(platemap$guide_id,
                           function(x) {
                             ifelse(x %in% mappings$new.NCR.id,
                                    mappings$antitag_label[mappings$new.NCR.id==x],
                                    ifelse(x %in% mappings$NCR.id,
                                           mappings$antitag_label[mappings$NCR.id==x],
                                           NA))
                           })
platemap$complementary <- complementary_nt[platemap$tag] == platemap$antitag

plate_samples <- unique(platemap[, c("guide_id", "conc", "group", "tag", "antitag", "complementary")])
plate_samples <- subset(plate_samples, guide_id != "NoProtein")
plate_samples$group <- with(plate_samples, paste0(group, " (", antitag, ")"))

# load plate reader data --------------------------------------------------

# read in file
plate_data <- readxl::read_xlsx(file.path(data_dir, data_fname),
                                range="A55:NW115")
plate_data$time <- plate_data$`Time [s]` / 60

# parse data by guide
plate_data <- lapply(seq(nrow(platemap)),
                     function(x) {
                       data.frame(time=plate_data$time,
                                  RFU=plate_data[[platemap$plate_cell[x]]],
                                  plate_cell=platemap$plate_cell[x],
                                  guide_id=platemap$guide_id[x],
                                  conc=platemap$conc[x],
                                  group=platemap$group[x])
                     })
plate_data <- do.call(rbind, plate_data)

# compare rates for 167fM -------------------------------------------------

guide_rates <- lapply(unique(mappings$NCR.id),
                      function(old_guide) {
                        new_guides <- mappings$new.NCR.id[mappings$NCR.id==old_guide]
                        tmp_data <- within(subset(plate_data,
                                                  guide_id %in% c(old_guide, new_guides) & 
                                                    conc=="167fM" & time >= 10),
                                           guide_id <- factor(guide_id,
                                                              levels=c(old_guide, new_guides)))
                        guide_rates <- lapply(levels(tmp_data$guide_id),
                                              function(x) {
                                                tmp_model <- lm(RFU ~ time,
                                                                data=subset(tmp_data,
                                                                            guide_id==x))
                                                tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                tmp_coef$guide_id <- x
                                                tmp_coef$term <- rownames(tmp_coef)
                                                return(subset(tmp_coef, term=="time"))
                                              })
                        guide_rates <- do.call(rbind, guide_rates)
                        guide_rates$p_diff <- sapply(guide_rates$guide_id,
                                                     function(x) {
                                                       if(x == old_guide) {
                                                         return(NA)
                                                       } else {
                                                         tmp_model <- lm(RFU ~ time * guide_id,
                                                                         data=subset(tmp_data,
                                                                                     guide_id %in% c(x, old_guide)))
                                                         tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                         return(tmp_coef$Pr...t..[grepl("time:", rownames(tmp_coef))])
                                                       }
                                                     })
                        return(guide_rates)
                      })
guide_rates <- do.call(rbind, guide_rates)
guide_rates <- dplyr::left_join(guide_rates, 
                                subset(plate_samples, conc=="167fM"), 
                                by="guide_id")
guide_rates$group <- factor(guide_rates$group,
                            levels=c("NCR_518 (A)", "NCR_564 (C)",
                                     "NCR_504 (G)", "NCR_579 (G)",
                                     "NCR_1362 (G)", "NCR_1379 (G)"))
guide_rates$tag <- factor(guide_rates$tag, levels=c("A", "G", "U", "C"))
guide_rates$guide_id <- factor(guide_rates$guide_id,
                               levels=unlist(lapply(levels(guide_rates$tag),
                                                    function(x) {
                                                      guide_rates$guide_id[guide_rates$tag==x]
                                                    })))

guide_rate_plot <- ggplot(guide_rates,
                          aes(x=guide_id, y=Estimate, 
                              ymin=Estimate + qnorm(0.025)*Std..Error,
                              ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col(aes(fill=tag)) + geom_errorbar(width=0.5, aes(col=complementary)) + 
  theme_classic() + facet_grid(~group, scales="free_x") + 
  xlab("") + ylab("RFU/min") + 
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5)) + 
  scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) + 
  labs(col="complementary\ntag") + 
  ggtitle("Fluorescence rates for 167 fM")


# compare rates for 0fM ---------------------------------------------------

bkgd_rates <- lapply(unique(mappings$NCR.id),
                     function(old_guide) {
                       new_guides <- mappings$new.NCR.id[mappings$NCR.id==old_guide]
                       tmp_data <- within(subset(plate_data,
                                                 guide_id %in% c(old_guide, new_guides) & 
                                                   conc=="0fM" & time >= 10),
                                          guide_id <- factor(guide_id,
                                                             levels=c(old_guide, new_guides)))
                       bkgd_rates <- lapply(levels(tmp_data$guide_id),
                                            function(x) {
                                              tmp_model <- lm(RFU ~ time,
                                                              data=subset(tmp_data,
                                                                          guide_id==x))
                                              tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                              tmp_coef$guide_id <- x
                                              tmp_coef$term <- rownames(tmp_coef)
                                              return(subset(tmp_coef, term=="time"))
                                            })
                       bkgd_rates <- do.call(rbind, bkgd_rates)
                       bkgd_rates$p_diff <- sapply(bkgd_rates$guide_id,
                                                   function(x) {
                                                     if(x == old_guide) {
                                                       return(NA)
                                                     } else {
                                                       tmp_model <- lm(RFU ~ time * guide_id,
                                                                       data=subset(tmp_data,
                                                                                   guide_id %in% c(x, old_guide)))
                                                       tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                       return(tmp_coef$Pr...t..[grepl("time:", rownames(tmp_coef))])
                                                     }
                                                   })
                       return(bkgd_rates)
                     })
bkgd_rates <- do.call(rbind, bkgd_rates)
bkgd_rates <- dplyr::left_join(bkgd_rates, 
                               subset(plate_samples, conc=="167fM"), 
                               by="guide_id")
bkgd_rates$group <- factor(bkgd_rates$group,
                           levels=c("NCR_518 (A)", "NCR_564 (C)",
                                    "NCR_504 (G)", "NCR_579 (G)",
                                    "NCR_1362 (G)", "NCR_1379 (G)"))
bkgd_rates$tag <- factor(bkgd_rates$tag, levels=c("A", "G", "U", "C"))
bkgd_rates$guide_id <- factor(bkgd_rates$guide_id,
                              levels=unlist(lapply(levels(bkgd_rates$tag),
                                                   function(x) {
                                                     bkgd_rates$guide_id[bkgd_rates$tag==x]
                                                   })))

bkgd_rate_plot <- ggplot(bkgd_rates,
                         aes(x=guide_id, y=Estimate, 
                             ymin=Estimate + qnorm(0.025)*Std..Error,
                             ymax=Estimate + qnorm(0.975)*Std..Error)) + 
  geom_col(aes(fill=tag)) + geom_errorbar(width=0.5, aes(col=complementary)) + 
  theme_classic() + facet_grid(~group, scales="free_x") + 
  xlab("") + ylab("RFU/min") + 
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5)) + 
  scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) + 
  labs(col="complementary\ntag") + 
  ggtitle("Fluorescence rates for 0 fM")
