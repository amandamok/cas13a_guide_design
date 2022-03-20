rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load computed slopes ----------------------------------------------------

guide_rate <- read.csv(file.path(data_dir, "guide_rate.csv"))
guide_rate$crRNA_spacer_basepairs <- as.factor(guide_rate$crRNA_spacer_basepairs)

noActivator_rate <- read.csv(file.path(data_dir, "noActivator_rate.csv"))
noActivator_rate$crRNA_spacer_basepairs <- as.factor(noActivator_rate$crRNA_spacer_basepairs)

# compute pvalues ---------------------------------------------------------

spacer_structure_antitag <- lapply(unique(guide_rate$antitag_pos1),
                                          function(x) {
                                            tmp_data <- subset(guide_rate,
                                                               crRNA_spacer_basepairs %in% c(0, 4, 10, 16) & 
                                                                 antitag_pos1 == x)
                                            tmp_model <- lm(Estimate ~ as.numeric(as.character(crRNA_spacer_basepairs)),
                                                            tmp_data)
                                            return(data.frame(summary(tmp_model)$coefficients,
                                                              antitag_pos1=x))
                                          })
spacer_structure_antitag <- do.call(rbind, spacer_structure_antitag)
spacer_structure_antitag_text <- data.frame(x=8, 
                                            y=rep(c(80, 70), each=4),
                                            text=c(paste("slope =",
                                                         round(spacer_structure_antitag$Estimate[2*(1:4)],
                                                               digits=3)),
                                                   paste("p =",
                                                         signif(spacer_structure_antitag$Pr...t..[2*(1:4)],
                                                                digits=3))),
                                            antitag_pos1=rep(unique(spacer_structure_antitag$antitag_pos1), 
                                                             times=2))

# generate plots ----------------------------------------------------------

suppl_figure_4A_noActivator <- ggplot(subset(noActivator_rate, 
                                             crRNA_spacer_basepairs %in% c(0, 4, 10, 16)),
                                      aes(x=as.numeric(as.character(crRNA_spacer_basepairs)), 
                                          y=Estimate)) +
  geom_violin(aes(fill=crRNA_spacer_basepairs),
              draw_quantiles=c(0.25, 0.75), linetype="dashed") +
  geom_violin(aes(fill=crRNA_spacer_basepairs), alpha=0.25, draw_quantiles=0.5) +
  geom_dotplot(aes(group=crRNA_spacer_basepairs), fill="grey35",
               binaxis="y", stackdir="center", binwidth=1) +
  geom_smooth(method=lm, formula=y~x) +
  theme_classic(base_size=10) + guides(fill="none") + 
  scale_x_continuous(breaks=c(0, 4, 10, 16)) +
  xlab("# basepaired positions in spacer") + ylab("RFU/min") + 
  ggtitle("activator-independent rate")


suppl_figure_4A_100fM <- ggplot(subset(guide_rate, 
                                       crRNA_spacer_basepairs %in% c(0, 4, 10, 16)),
                                aes(x=as.numeric(as.character(crRNA_spacer_basepairs)), 
                                    y=Estimate)) +
  geom_violin(aes(fill=crRNA_spacer_basepairs),
              draw_quantiles=c(0.25, 0.75), linetype="dashed") +
  geom_violin(aes(fill=crRNA_spacer_basepairs), alpha=0.25, draw_quantiles=0.5) +
  geom_dotplot(aes(group=crRNA_spacer_basepairs), fill="grey35",
               binaxis="y", stackdir="center", binwidth=1) +
  geom_smooth(method=lm, formula=y~x) +
  theme_classic(base_size=10) + guides(fill="none") + 
  scale_x_continuous(breaks=c(0, 4, 10, 16)) +
  xlab("# basepaired positions in spacer") + ylab("RFU/min") + 
  ggtitle("activator-dependent rate") + facet_grid(~antitag_pos1) + 
  geom_text(data=spacer_structure_antitag_text, 
            aes(x=x, y=y, label=text), hjust=0.5, size=2.5)

suppl_figure_4A <- suppl_figure_4A_noActivator / suppl_figure_4A_100fM

ggsave(filename=file.path(figure_dir, "suppl_figure_4A.pdf"),
       plot=suppl_figure_4A,
       device="pdf", width=6.5, height=4, units="in")
