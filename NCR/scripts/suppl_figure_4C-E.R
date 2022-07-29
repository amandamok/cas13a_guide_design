rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))

# perform correlation tests -----------------------------------------------

# GC content
GC_content_test <- with(guide_rate,
                        cor.test(GC_content, Estimate, method="spearman"))

# gRNA MFE
gRNA_MFE_test <- with(guide_rate,
                      cor.test(gRNA_MFE, Estimate, method="spearman"))

# U context
U_context_data <- data.frame(Estimate=rep(guide_rate$Estimate, times=2),
                             U_context=c(guide_rate$U_context, 
                                         guide_rate$U_context_unstructured),
                             structure=rep(c("structured + unstructured", 
                                             "unstructured"),
                                           each=nrow(guide_rate)))
U_context_test <- list(both=cor.test(guide_rate$Estimate, 
                                    guide_rate$U_context),
                      unstructured=cor.test(guide_rate$Estimate, 
                                            guide_rate$U_context_unstructured))

# generate plots ----------------------------------------------------------

figure_S4C <- ggplot(guide_rate, aes(x=GC_content, y=Estimate)) + 
  geom_point(col="grey35") + theme_classic(base_size=8) +
  geom_smooth(method="lm", formula=y~x) + 
  # ggtitle(bquote(rho == .(signif(GC_content_test$estimate, 3)))) + 
  xlab("spacer %GC") + ylab("activator-dependent activity\nRFU/min")

figure_S4D <- ggplot(guide_rate, aes(x=gRNA_MFE, y=Estimate)) + 
  geom_point(col="grey35") + theme_classic(base_size=8) +
  geom_smooth(method="lm", formula=y~x) + 
  # ggtitle(bquote(rho == .(signif(gRNA_MFE_test$estimate, 3)))) + 
  xlab("crRNA MFE\nkcal/mol") + ylab("activator-dependent activity\nRFU/min")

figure_S4E <- ggplot(U_context_data, aes(x=U_context, y=Estimate)) + 
  geom_point() + geom_smooth(method=lm, formula=y~x) + 
  facet_grid(~structure, scales="free_x") + 
  theme_classic(base_size=8) + 
  xlab("%U in target neighborhood") + 
  ylab("activator-dependent rate\n(RFU/min)")

ggsave(filename=file.path(figure_dir, "suppl_figure_4C.pdf"),
       plot=figure_S4C,
       device="pdf", width=3, height=2, units="in")
ggsave(filename=file.path(figure_dir, "suppl_figure_4D.pdf"),
       plot=figure_S4D,
       device="pdf", width=3, height=2, units="in")
ggsave(filename=file.path(figure_dir, "suppl_figure_4E.pdf"),
       plot=figure_S4E,
       device="pdf", width=6.5, height=2.25, units="in")
