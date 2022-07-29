rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data", "supplementary_data")
figure_dir <- file.path(project_dir, "figures")

# load platemap -----------------------------------------------------------

platemap_fname <- file.path(data_dir, "Human RNA Screen Platemap.xlsx")
platemap <- readxl::read_xlsx(platemap_fname,
                              range="B2:Y18")

platemap <- data.frame(plate_row=rep(LETTERS[seq(nrow(platemap))],
                                     times=ncol(platemap)),
                       plate_col=rep(seq(ncol(platemap)), each=nrow(platemap)),
                       guide=unlist(platemap))
platemap$plate_well <- with(platemap, paste0(plate_row, plate_col))
platemap <- subset(platemap, !is.na(platemap$guide))
platemap$activator <- ifelse(platemap$plate_row %in% LETTERS[1:8],
                             "+", "-")
platemap$sample <- with(platemap, paste(guide, activator, sep="_"))

# load plate reader data --------------------------------------------------

plate_data_fname <- file.path(data_dir, "HumanRNA_1.xlsx")
plate_data <- readxl::read_xlsx(plate_data_fname, range="A55:NW115")
plate_data$time <- round(plate_data$`Time [s]`/60)

plate_data <- lapply(platemap$plate_well,
                     function(x) {
                       data.frame(time=plate_data$time,
                                  RFU=plate_data[[x]],
                                  plate_well=x)
                     })
plate_data <- do.call(rbind, plate_data)

plate_data <- plyr::join(plate_data, 
                         platemap[, c("plate_well", "guide", "activator", "sample")],
                         by="plate_well")

# compute slopes ----------------------------------------------------------

human_rna_slopes <- lapply(unique(plate_data$guide),
                           function(x) {
                             tmp_data <- subset(plate_data, guide==x)
                             tmp_model <- lm(RFU ~ time * activator, tmp_data)
                             tmp_coef <- data.frame(summary(tmp_model)$coefficients,
                                                    guide=x)
                             tmp_coef$coef <- rownames(tmp_coef)
                             return(tmp_coef)
                           })
human_rna_slopes <- do.call(rbind, human_rna_slopes)
human_rna_slopes <- subset(human_rna_slopes, grepl("time:", human_rna_slopes$coef))

human_rna_slopes$type <- "expt"
human_rna_slopes$type[human_rna_slopes$guide %in% c("No Protein", "Apo Protein", "612")] <- "control"
human_rna_slopes$text <- ""
human_rna_slopes$text[p.adjust(human_rna_slopes$Pr...t.., method="fdr")<0.05] <- "*"
human_rna_slopes$guide[human_rna_slopes$guide=="Apo Protein"] <- "RNP-only"
human_rna_slopes$guide[human_rna_slopes$guide=="612"] <- "crRNA control"

# generate plot -----------------------------------------------------------

figure_S8B <- ggplot(human_rna_slopes, 
       aes(x=factor(guide, levels=c("No Protein", "RNP-only", "crRNA control",
                                    sort(as.numeric(guide)))), 
           y=Estimate, fill=type)) + 
  geom_col() + theme_classic(base_size=8) + scale_fill_manual(values=c("grey", "black")) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), 
        legend.position="none") + 
  geom_errorbar(aes(ymin=Estimate+qnorm(0.025)*Std..Error, 
                    ymax=Estimate-qnorm(0.025)*Std..Error)) + 
  geom_text(aes(y=Estimate - qnorm(0.025) * Std..Error + 0.005, label=text, col="red")) + 
  xlab("") + ylab("activator-dependent rate\n(RFU/min)")

ggsave(file=file.path(figure_dir, "suppl_figure_8B.pdf"),
       plot=figure_S8B,
       device="pdf", width=6.5, height=2.5, units="in")
