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

concentrations <- c("no_protein", "no_activator",
                    "125", "62", "31", "16", "8", "4", "2", "1", "0")
platemap <- data.frame(plate_row=rep(LETTERS[1:4], each=22),
                       plate_col=rep(c(2:23), times=4),
                       conc=rep(rep(concentrations, each=2), times=4))
platemap$multiplex <- ifelse(platemap$plate_row %in% LETTERS[1:2],
                             "8g", "32g")
platemap$well <- with(platemap, paste(plate_row, plate_col, sep=""))

data_fname <- file.path(data_dir, "EM_SecondaryLod_Rep4_06242022_75_55g.xlsx")
raw_data <- readxl::read_xlsx(data_fname, range="A76:NW136")
raw_data <- reshape2::melt(raw_data, id.vars=1:3, variable.name="well")
raw_data$time <- raw_data$`Time [s]`/60
raw_data$value <- as.numeric(raw_data$value)
raw_data <- subset(raw_data, well %in% platemap$well)
raw_data <- droplevels(raw_data)
raw_data <- plyr::join(raw_data, platemap, by="well")

endpoint_time <- max(raw_data$time)

# normalize data ----------------------------------------------------------

normalized_data <- lapply(split(raw_data, raw_data$well),
                          function(x) {
                            which_row <- which(round(x$time) == 6)
                            baseline <- x$value[which_row]
                            return(within(x, {
                              RFU_norm <- value/baseline
                            }))
                          })
normalized_data <- do.call(rbind, normalized_data)

# exploratory analyses ----------------------------------------------------

ggplot(normalized_data,
       aes(x=time, y=value, group=well,
           col=factor(conc, levels=unique(platemap$conc)))) + 
  geom_line(alpha=0.5) + theme_classic() + facet_grid(~multiplex) + 
  xlab("time (min)") + ylab("RFU") + labs(col="") 

ggplot(normalized_data,
       aes(x=time, y=RFU_norm, group=well,
           col=factor(conc, levels=unique(platemap$conc)))) + 
  geom_line(alpha=0.5) + theme_classic() + facet_grid(~multiplex) + 
  xlab("time (min)") + ylab(TeX("$RFU_t / RFU_{t=6}$")) + labs(col="")

ggplot(subset(normalized_data, time==endpoint_time),
       aes(x=well, y=RFU_norm, fill=plate_row)) + 
  geom_col() + theme_classic() + 
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) + 
  xlab("") + ylab(TeX("$RFU_t / RFU_{t=6}$")) + guides(fill="none")


# test endpoint fluorescence ----------------------------------------------

endpoint_data <- subset(normalized_data, time==max(normalized_data$time))
endpoint_data$conc <- factor(endpoint_data$conc, levels=concentrations)

endpoint_pvalues <- endpoint_data %>% 
  group_by(multiplex) %>% 
  t_test(RFU_norm ~ conc, ref.group="no_activator", alternative="less", p.adjust.method="none") %>%
  add_xy_position(x="conc") %>%
  add_significance()

# generate plot -----------------------------------------------------------

endpoint_means <- within(aggregate(RFU_norm ~ conc + multiplex, 
                                   endpoint_data, FUN=mean),
                         std_dev <- aggregate(RFU_norm ~ conc + multiplex, 
                                              endpoint_data, FUN=sd)$RFU_norm)
endpoint_means$pvalue <- endpoint_pvalues$p[prodlim::row.match(endpoint_means[, c("conc", "multiplex")],
                                                               endpoint_pvalues[, c("group2", "multiplex")])]
endpoint_means$pvalue_label <- ifelse(endpoint_means$pvalue < 0.05,
                                      "p < 0.05", "ns")
endpoint_means$pvalue_label[is.na(endpoint_means$pvalue_label)] <- "ns"

ggplot(endpoint_means, 
       aes(x=factor(conc, levels=concentrations), y=RFU_norm)) + 
  geom_col(aes(fill=pvalue_label)) + theme_classic(base_size=10) + 
  geom_errorbar(aes(ymin=RFU_norm + qnorm(0.025)*std_dev,
                    ymax=RFU_norm + qnorm(0.975)*std_dev), width=0.5) + 
  facet_grid(~factor(multiplex, levels=c("8g", "32g"))) + 
  geom_point(data=endpoint_data, shape=1) + 
  stat_pvalue_manual(endpoint_pvalues) + labs(fill="") + 
  scale_fill_manual(values=c("p < 0.05"="red", "ns"="grey35")) + 
  xlab(TeX("SARS-CoV-2 RNA (copies/$\\mu$L)")) + 
  ylab(TeX("Fluorescence ($RFU/RFU_6$)")) + 
  theme(panel.spacing.x=unit(3, "lines")) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
