rm(list=ls())

library(ggplot2)
library(ggpattern)

plate_nrow <- 16
plate_ncol <- 24

complementary_nt <- c("A", "U", "C", "G")
names(complementary_nt) <- c("U", "A", "G", "C")

time_delay <- 10 # timepoint to start calculating slopes

# filepaths
guide_map_fname <- "../NCR inventory - test_alternative_tag.csv"
plate_map_fname <- "Dylan - Platemap 011422 Anti-Tag gBlocks Screen.xlsx - Sheet1.csv"
raw_data_fname <- "Anti-tag 504 564 @ 250 fM Cas13 LoD (Modified)_20220114.xlsx - Result sheet.csv"
primary_screen_slopes_fname <- "../../guide_rate.csv"

# load guide map ----------------------------------------------------------

guide_map <- read.csv(guide_map_fname)

# load plate map ----------------------------------------------------------

plate_map <- read.csv(plate_map_fname)
plate_map <- plate_map[-1, -1]
plate_map <- data.frame(plate_well = paste0(rep(LETTERS[seq(plate_nrow)], times=plate_ncol),
                                            rep(seq(plate_ncol), each=plate_nrow)),
                        sample=unlist(plate_map),
                        stringsAsFactors=F)
# account for merged cells
for(x in seq(plate_ncol)[(seq(plate_ncol) %% 2) == 0]) { # even columns
  if(x %in% c(2:12)) {
    match_rows <- 1:8
  } else {
    if(x %in% 13:22) {
      match_rows <- 1:6
    } else {
      match_rows <- 3:6
    }
  }
  old_cells <- paste0(LETTERS[match_rows], x-1)
  new_cells <- paste0(LETTERS[match_rows], x)
  plate_map$sample[match(new_cells, plate_map$plate_well)] <- plate_map$sample[match(old_cells, plate_map$plate_well)]
}
plate_map <- subset(plate_map, !is.na(plate_map$sample) & !(plate_map$sample == ""))

# guide features ----------------------------------------------------------

# add guides that weren't screened
missing_guides <- sapply(unique(guide_map$NCR.id),
                         function(x) {
                           ifelse(!any(grepl(x, plate_map$sample)),
                                  x, NA)
                         })
missing_guides <- missing_guides[!is.na(missing_guides)]
missing_samples <- paste(rep(missing_guides, each=2),
                         rep(c("No activator", "100 fM vRNA"), times=length(missing_guides)))

# initialize data.frame
sample_features <- data.frame(sample=c(unique(plate_map$sample),
                                       missing_samples))

# label activator 
sample_features$activator <- sapply(sample_features$sample,
                                    function(x) {
                                      ifelse(x == "No Protein",
                                             "No protein",
                                             ifelse(grepl("fM", x),
                                                    "+ activator",
                                                    "- activator"))
                                    })
sample_features$activator_type <- sapply(sample_features$sample,
                                         function(x) {
                                           ifelse(x == "No Protein",
                                                  "No protein",
                                                  ifelse(grepl("fM", x),
                                                         ifelse(grepl("gBlocks", x),
                                                                "250 fM gBlocks",
                                                                "100 fM vRNA"),
                                                         "No activator"))
                                         })

# label guide
sample_features$guide_id <- sapply(sample_features$sample,
                                   function(x) {
                                     if(x == "No Protein") {
                                       return("No Protein") 
                                     } else {
                                       return(strsplit(x, split=" ")[[1]][1])
                                     }
                                   })
sample_features$old_guide <- sample_features$guide_id %in% guide_map$NCR.id

# label (old) comparison guide
sample_features$group <- sapply(sample_features$guide_id,
                                function(x) {
                                  if(x %in% guide_map$new.NCR.id) {
                                    return(guide_map$NCR.id[match(x, guide_map$new.NCR.id)])
                                  } else {
                                    return(x)
                                  }
                                })

# label tag
sample_features$tag <- guide_map$new_tag_pos1[match(sample_features$guide_id, guide_map$new.NCR.id)]
sample_features$tag[sample_features$old_guide] <- "C"

# label anti-tag
sample_features$antitag <- guide_map$antitag_label[match(sample_features$group,
                                                         guide_map$NCR.id)]

# group label: original guide & anti-tag
sample_features$group_label <- paste0(sample_features$group, " (", sample_features$antitag, ")")
sample_features$group_label[sample_features$group_label == "NA (NA)"] <- "No Protein"

# label complementarity between tag & anti-tag
sample_features$complementary <- (sample_features$tag == complementary_nt[sample_features$antitag])

# load raw plate reader data ----------------------------------------------

# load data corresponding to 10 nm (excitation) / 20 nm (emission)
raw_data <- read.csv(raw_data_fname, header=T, skip=75, nrows=50)

# reshape data
plate_data <- lapply(seq(nrow(plate_map)),
                     function(x) {
                       plate_well <- plate_map$plate_well[x]
                       data.frame(time=raw_data$Time..s.,
                                  RFU=raw_data[[plate_well]],
                                  well=plate_well,
                                  sample=plate_map$sample[x])
                     })
plate_data <- do.call(rbind, plate_data)
plate_data$time <- plate_data$time/60 # convert time from seconds to minutes

# calculate slopes --------------------------------------------------------

sample_slopes <- lapply(unique(plate_data$sample),
                        function(tmp_sample) {
                          tmp_data <- subset(plate_data, plate_data$sample==tmp_sample)
                          tmp_model <- lm(RFU ~ time, 
                                          subset(tmp_data, time >= time_delay))
                          which_row <- which(rownames(summary(tmp_model)$coefficients) == "time")
                          return(summary(tmp_model)$coefficients[which_row,])
                        })
sample_slopes <- data.frame(sample=unique(plate_data$sample),
                            do.call(rbind, sample_slopes))
colnames(sample_slopes) <- c("sample", "Estimate", "Std.Error", "t.value", "p.value")

# add slopes from primary screen ------------------------------------------

primary_screen_slopes <- read.csv(primary_screen_slopes_fname)

# subset to needed guides
missing_guides_slopes <- subset(primary_screen_slopes, 
                                NCR.id %in% missing_guides,
                                select=c(NCR.id, Estimate, Std.Error, t.value, p.value))

# prepare for joining
colnames(missing_guides_slopes) <- c("sample", "Estimate", "Std.Error", "t.value", "p.value")
missing_guides_slopes$sample <- paste(missing_guides_slopes$sample, "100 fM vRNA")

# add slopes to sample_features -------------------------------------------

# merge slopes from this expt & primary screen
all_slopes <- rbind(sample_slopes, missing_guides_slopes)

# add sample features
all_slopes <- merge(x=all_slopes, 
                    y=as.data.frame(unclass(sample_features), stringsAsFactors=T), 
                    by="sample", all.y=T)

# plot results ------------------------------------------------------------

all_slopes$label_complementary <- ifelse(all_slopes$complementary,
                                         "***", "")
all_slopes$label_vRNA <- ifelse(all_slopes$activator_type == "100 fM vRNA",
                                "100 fM", "")
all_slopes$group_label <- factor(as.character(all_slopes$group_label),
                                 levels=levels(all_slopes$group_label)[c(6, 1, 2, 3, 4, 5, 7)])
all_slopes$tag <- factor(as.character(all_slopes$tag),
                         levels=c("A", "U", "G", "C"))

all_slopes$Estimate[is.na(all_slopes$Estimate)] <- 0
all_slopes$Std.Error[is.na(all_slopes$Std.Error)] <- 0

ggplot(all_slopes, 
       aes(x=tag, y=Estimate, 
           ymin=Estimate + qnorm(0.05)*Std.Error,
           ymax=Estimate + qnorm(0.95)*Std.Error,
           fill=tag)) + 
  theme_classic() + geom_col() + geom_errorbar() + 
  geom_text(aes(label=label_vRNA, y=Estimate + 2*Std.Error)) + 
  geom_text(aes(label=label_complementary, y=Estimate + 2.5*Std.Error)) + 
  # facet_wrap(group_label~activator, scales="free", nrow=2, dir="v", drop=T) + 
  facet_wrap(vars(group_label, activator), scales="free", drop=T, nrow=2, dir="v") + 
  geom_hline(yintercept=0) + xlab("") + ylab("RFU/min") + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
