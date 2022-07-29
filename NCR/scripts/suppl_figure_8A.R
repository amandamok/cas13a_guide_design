rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data", "supplementary_data")
figure_dir <- file.path(project_dir, "figures")

# load platemap -----------------------------------------------------------

platemap_fname <- file.path(data_dir, "LeaveOneOut Platemap.xlsx")
platemap <- readxl::read_xlsx(platemap_fname, range="B2:Y18")

platemap <- data.frame(plate_row=rep(LETTERS[seq(nrow(platemap))],
                                     times=ncol(platemap)),
                       plate_col=rep(seq(ncol(platemap)), each=nrow(platemap)),
                       sample=unlist(platemap))
platemap$plate_well <- with(platemap, paste0(plate_row, plate_col))
platemap <- subset(platemap, !is.na(platemap$sample))

pool_4 <- unique(grep("\\+", platemap$sample, value=T))
platemap$pool <- sapply(platemap$sample,
                        function(x) {
                          if(any(grepl(x, pool_4))) {
                            return(pool_4[grepl(x, pool_4)])
                          } else {
                            return(x)
                          }
                        })

# load plate reader data --------------------------------------------------

plate_data_fname <- file.path(data_dir, "LeaveOneOut.xlsx")
plate_data <- readxl::read_xlsx(plate_data_fname, range="A57:NW117")
plate_data$time <- round(plate_data$`Time [s]`/60)

plate_data <- lapply(platemap$plate_well,
                     function(x) {
                       data.frame(time=plate_data$time,
                                  RFU=plate_data[[x]],
                                  plate_well=x)
                     })
plate_data <- do.call(rbind, plate_data)

plate_data <- plyr::join(plate_data,
                         platemap[, c("plate_well", "sample", "pool")],
                         by="plate_well")

# compute slopes ----------------------------------------------------------

pool4_tests <- lapply(pool_4,
                      function(x) {
                        tmp_data <- subset(plate_data, pool==x & time >= 20)
                        tmp_data$sample <- relevel(as.factor(tmp_data$sample), ref=x)
                        tmp_model <- lm(RFU ~ time * sample, tmp_data)
                        tmp_anova <- anova(tmp_model)
                        tmp_slopes <- emmeans::emtrends(tmp_model, "sample", var="time")
                        tmp_pairs <- data.frame(pairs(tmp_slopes), pool=x)
                        tmp_pairs <- subset(tmp_pairs, grepl("\\(", tmp_pairs$contrast))
                        tmp_slopes <- data.frame(tmp_slopes)
                        tmp_slopes$p_contrast <- tmp_pairs$p.value[match(paste0("(", x, ") - ", 
                                                                                tmp_slopes$sample),
                                                                         tmp_pairs$contrast)]
                        return(tmp_slopes)
                      })
pool4_tests <- do.call(rbind, pool4_tests)

loo_control_tests <- data.frame(emmeans::emtrends(lm(RFU ~ time * sample, 
                                                     subset(plate_data, !(pool %in% pool_4) & 
                                                              time >= 20)), 
                                                  "sample", var="time"),
                                p_contrast = NA)

loo_plot_data <- rbind(pool4_tests, loo_control_tests)

# generate plot -----------------------------------------------------------

loo_plot_data$pool <- sapply(as.character(loo_plot_data$sample),
                             function(x) {
                               if(any(grepl(x, pool_4))) {
                                 return(pool_4[grepl(x, pool_4)])
                               } else {
                                 return(ifelse(x %in% pool_4, x, "control"))
                               }
                             })
levels(loo_plot_data$sample)[levels(loo_plot_data$sample)=="No Act"] <- "RNP-only"
levels(loo_plot_data$sample)[levels(loo_plot_data$sample)=="612"] <- "crRNA control"
loo_plot_data$sample <- sapply(as.character(loo_plot_data$sample),
                               function(x) {
                                 ifelse(grepl("\\+", x) | 
                                          (x %in% c("crRNA control", "8plex", "RNP-only", "No Protein")), 
                                        x, paste0("-", x))
                               })
loo_plot_data$sample <- factor(loo_plot_data$sample, levels=loo_plot_data$sample)
loo_plot_data$text <- ""
loo_plot_data$text[loo_plot_data$p_contrast < 0.05] <- "*"

figure_S8A <- ggplot(loo_plot_data, aes(x=sample, y=time.trend, fill=pool)) + geom_col() + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL)) + geom_hline(yintercept=0) + 
  geom_text(aes(y=time.trend+0.5, label=text)) + 
  theme_classic(base_size=8) + xlab("") + ylab("RFU/min") + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.position="none")

ggsave(file=file.path(figure_dir, "suppl_figure_8A.pdf"),
       plot=figure_S8A,
       device="pdf", width=6.5, height=2.5, units="in")
