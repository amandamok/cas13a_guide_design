rm(list=ls())

library(here)
library(ggplot2)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data", "supplementary_data")
figure_dir <- file.path(project_dir, "figures")

# generate platemap -------------------------------------------------------

variants_platemap <- data.frame(plate_row=rep(LETTERS[1:4], each=16),
                                plate_col=rep(1:16, times=4))
variants_platemap$plate_well <- with(variants_platemap, paste0(plate_row, plate_col))
variants_platemap <- subset(variants_platemap,
                            !((plate_row %in% LETTERS[c(2, 4)]) & (plate_col %in% (2*seq(8)-1))))
variants_platemap$multiplex <- ifelse(variants_platemap$plate_row %in% LETTERS[1:2],
                                      8, 32)
variants_platemap$multiplex <- as.factor(variants_platemap$multiplex)
variant_samples <- c("No\nprotein", "No\nactivator", "WT", "Alpha", 
                     "Beta", "Gamma", "Kappa", "Delta")
variants_platemap$variant <- ""
for(x in seq_along(variant_samples)) {
  variants_platemap$variant[variants_platemap$plate_col %in% c(2*x-1, 2*x)] <- variant_samples[x]
}
variants_platemap$variant <- factor(variants_platemap$variant, levels=variant_samples)

# load data ---------------------------------------------------------------

variants_data <- readxl::read_xlsx(file.path(data_dir, "variant_detection_EM_28092021_data.xlsx"),
                                     range="A84:NW144")
variants_data$time <- round(variants_data$`Time [s]`/60)

variants_data <- lapply(variants_platemap$plate_well,
                        function(x) {
                          data.frame(time=variants_data$time,
                                     RFU=as.numeric(variants_data[[x]]),
                                     plate_well=x)
                        })
variants_data <- do.call(rbind, variants_data)

variants_data <- plyr::join(variants_data, 
                            variants_platemap[, c("plate_well", "multiplex", "variant")],
                            by="plate_well")

# calculate slopes --------------------------------------------------------

variant_slopes <- lapply(variant_samples,
                         function(x) {
                           tmp_slopes <- lapply(c(8, 32),
                                                function(y) {
                                                  tmp_data <- subset(variants_data,
                                                                     variant==x & multiplex==y & time >= 10)
                                                  tmp_model <- lm(RFU ~ time, tmp_data)
                                                  tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                                  tmp_coef$coef <- rownames(tmp_coef)
                                                  tmp_coef$variant <- x
                                                  tmp_coef$multiplex <- y
                                                  return(tmp_coef)
                                                })
                           tmp_slopes <- do.call(rbind, tmp_slopes)
                           return(tmp_slopes)
                         })
variant_slopes <- do.call(rbind, variant_slopes)
variant_slopes <- subset(variant_slopes, coef=="time")
variant_slopes$variant <- factor(variant_slopes$variant, levels=variant_samples)
variant_slopes$multiplex <- factor(variant_slopes$multiplex, levels=c(8, 32))

# compare slopes between 8-plex and 32-plex -------------------------------

variant_compare_slopes <- lapply(variant_samples,
                                 function(x) {
                                   tmp_data <- subset(variants_data,
                                                      variant==x & time >= 10)
                                   tmp_model <- lm(RFU ~ time*multiplex, tmp_data)
                                   tmp_coef <- data.frame(summary(tmp_model)$coefficients)
                                   tmp_coef$coef <- rownames(tmp_coef)
                                   tmp_coef$variant <- x
                                   return(tmp_coef)
                                 })
variant_compare_slopes <- do.call(rbind, variant_compare_slopes)
variant_compare_slopes <- subset(variant_compare_slopes, 
                                 variant_compare_slopes$coef=="time:multiplex32")

# add labels for ggprism::add_pvalue()
variant_compare_slopes$diff <- cut(variant_compare_slopes$Pr...t..,
                                   breaks=c(0, 0.001, 0.01, 0.05, 1),
                                   labels=c("***", "**", "*", "ns"))
variant_compare_slopes$group1 <- 8
variant_compare_slopes$group2 <- 32
variant_compare_slopes$y.position <- sapply(seq(nrow(variant_compare_slopes)),
                                            function(x) {
                                              max(subset(variant_slopes, 
                                                         variant==variant_compare_slopes$variant[x])$Estimate)
                                            })
variant_compare_slopes$y.position <- variant_compare_slopes$y.position + 50
variant_compare_slopes$x <- as.numeric(factor(variant_compare_slopes$variant,
                                              levels=variant_samples))
variant_compare_slopes$xmin <- variant_compare_slopes$x - 0.2
variant_compare_slopes$xmax <- variant_compare_slopes$x + 0.2

# generate plot -----------------------------------------------------------

figure_6D <- ggplot(variant_slopes, aes(x=variant, y=Estimate, 
                                        ymin=Estimate+qnorm(0.05)*Std..Error,
                                        ymax=Estimate+qnorm(0.95)*Std..Error)) + 
  geom_col(aes(fill=multiplex), position="dodge", width=0.5) + 
  geom_errorbar(aes(group=multiplex), position=position_dodge(), width=0.5) + 
  theme_classic(base_size=8) + xlab("") + ylab("RFU/min") + 
  scale_fill_manual(values=c("8"="grey35", "32"="red"),
                    labels=c("8", "32")) + 
  labs(fill="multiplex\npool")

ggsave(filename=file.path(figure_dir, "figure_6D.pdf"),
       plot=figure_6D,
       device="pdf", width=6.5, height=2, units="in")
