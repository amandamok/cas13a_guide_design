rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load platemap -----------------------------------------------------------

loo_largepool_platemap <- readxl::read_xlsx(file.path(project_dir, "data", 
                                                      "supplementary_data", 
                                                      "LeaveOneOut_LargePool Platemap.xlsx"),
                                            sheet="Variants Platemap",
                                            range="B2:O10")
loo_largepool_platemap <- data.frame(plate_row=rep(LETTERS[seq(nrow(loo_largepool_platemap))],
                                                   times=ncol(loo_largepool_platemap)),
                                     plate_col=rep(seq(ncol(loo_largepool_platemap)),
                                                   each=nrow(loo_largepool_platemap)),
                                     multiplex=unlist(loo_largepool_platemap))
loo_largepool_platemap$well <- with(loo_largepool_platemap, paste0(plate_row, plate_col))
loo_largepool_platemap <- subset(loo_largepool_platemap, multiplex != "<NA>")
loo_largepool_conc <- c("No Protein", "RNP-only", "30 fM", "10 fM", "3 fM", "1 fM", "300 aM")
loo_largepool_platemap$conc <- ""
for(x in seq_along(loo_largepool_conc)) {
  tmp_col <- c(2*x-1, 2*x)
  tmp_conc <- loo_largepool_conc[x]
  loo_largepool_platemap$conc[loo_largepool_platemap$plate_col %in% tmp_col] <- tmp_conc
}
loo_largepool_platemap$sample <- with(loo_largepool_platemap, 
                                      paste(multiplex, conc, sep="_"))

# load data ---------------------------------------------------------------

# use findit_75 for gain=75
loo_largepool_data <- readxl::read_xlsx(file.path(project_dir, "data",
                                                  "supplementary_data",
                                                  "LeaveOneOut_LargePool.xlsx"),
                                        sheet="Result sheet", range="B101:GE161")
loo_largepool_data <- lapply(seq(nrow(loo_largepool_platemap)),
                             function(x) {
                               data.frame(time=loo_largepool_data$`Time [s]`,
                                          RFU=loo_largepool_data[[loo_largepool_platemap$well[x]]],
                                          well=loo_largepool_platemap$well[x],
                                          multiplex=loo_largepool_platemap$multiplex[x],
                                          conc=loo_largepool_platemap$conc[x],
                                          sample=loo_largepool_platemap$sample[x])
                             })
loo_largepool_data <- do.call(rbind, loo_largepool_data)
loo_largepool_data$time <- loo_largepool_data$time / 60

# compute slopes ----------------------------------------------------------

loo_largepool_rates <- lapply(unique(loo_largepool_data$sample),
                              function(x) {
                                tmp_data <- subset(loo_largepool_data, sample==x)
                                tmp_model <- lm(RFU ~ time, tmp_data)
                                tmp_coef <- data.frame(summary(tmp_model)$coefficients,
                                                       sample=x)
                                tmp_coef$coef <- rownames(tmp_coef)
                                return(tmp_coef)
                              })
loo_largepool_rates <- do.call(rbind, loo_largepool_rates)
loo_largepool_rates <- subset(loo_largepool_rates, coef=="time")
loo_largepool_rates$conc <- matrix(unlist(strsplit(loo_largepool_rates$sample, split="_")), 
                                   byrow=T, ncol=2)[,2]
loo_largepool_rates$multiplex <- matrix(unlist(strsplit(loo_largepool_rates$sample, split="_")),
                                        byrow=T, ncol=2)[,1]

# compute p-values --------------------------------------------------------

loo_largepool_detectable <- lapply(unique(loo_largepool_platemap$multiplex),
                                   function(tmp_pool) {
                                     if(tmp_pool == "No Protein") { return(NULL) }
                                     do.call(rbind,
                                             lapply(loo_largepool_conc[-c(1:2)],
                                                    function(tmp_conc) {
                                                      tmp_data <- subset(loo_largepool_data,
                                                                         multiplex==tmp_pool &
                                                                           conc %in% c("RNP-only", tmp_conc))
                                                      tmp_data$conc <- relevel(as.factor(tmp_data$conc),
                                                                               ref="RNP-only")
                                                      tmp_model <- lm(RFU ~ time * conc, tmp_data)
                                                      tmp_slopes <- data.frame(summary(tmp_model)$coefficients,
                                                                               multiplex=tmp_pool,
                                                                               conc=tmp_conc)
                                                      tmp_slopes$coef <- rownames(tmp_slopes)
                                                      return(tmp_slopes)
                                                    }))
                                   })
loo_largepool_detectable <- do.call(rbind, loo_largepool_detectable)
loo_largepool_detectable <- subset(loo_largepool_detectable, 
                                   grepl("time:", loo_largepool_detectable$coef))
loo_largepool_detectable$label <- ifelse(loo_largepool_detectable$Pr...t.. < 0.05,
                                         "*", "")
## all multiplex pools are detectable above RNP-only at all concentrations
loo_largepool_mismatch <- lapply(loo_largepool_conc[-(1:2)],
                                 function(tmp_conc) {
                                   do.call(rbind,
                                           lapply(c(8, 32),
                                                  function(tmp_pool) {
                                                    if(tmp_pool == 8) {
                                                      tmp_pools <- c("8 Guides", "7 Guides")
                                                    } else {
                                                      tmp_pools <- c("32 Guides", "31 Guides")
                                                    }
                                                    tmp_data <- subset(loo_largepool_data,
                                                                       conc==tmp_conc & 
                                                                         multiplex %in% tmp_pools)
                                                    tmp_data$multiplex <- relevel(as.factor(tmp_data$multiplex),
                                                                                  ref=tmp_pools[1])
                                                    tmp_model <- lm(RFU ~ time*multiplex, tmp_data)
                                                    tmp_coef <- data.frame(summary(tmp_model)$coefficients,
                                                                           multiplex_label = tmp_pools[1],
                                                                           conc=tmp_conc)
                                                    tmp_coef$coef <- rownames(tmp_coef)
                                                    return(tmp_coef)
                                                  }))
                                 })
loo_largepool_mismatch <- do.call(rbind, loo_largepool_mismatch)
loo_largepool_mismatch <- subset(loo_largepool_mismatch,
                                 grepl("time:", loo_largepool_mismatch$coef))
loo_largepool_mismatch$text <- ifelse(loo_largepool_mismatch$Pr...t.. < 0.05,
                                      "*", "")
loo_largepool_mismatch$text_height <- sapply(seq(nrow(loo_largepool_mismatch)),
                                             function(x) {
                                               tmp_pool <- loo_largepool_mismatch$multiplex[x]
                                               tmp_conc <- loo_largepool_mismatch$conc[x]
                                               if(tmp_pool == "8 Guides") {
                                                 tmp_pools <- c("8 Guides", "7 Guides")
                                               } else {
                                                 tmp_pools <- c("32 Guides", "31 Guides")
                                               }
                                               tmp_data <- subset(loo_largepool_rates,
                                                                  multiplex %in% tmp_pools &
                                                                    conc == tmp_conc)
                                               return(max(tmp_data$Estimate))
                                             })
loo_largepool_mismatch$multiplex_label <- relevel(as.factor(loo_largepool_mismatch$multiplex_label),
                                                  ref="8 Guides")
loo_largepool_mismatch <- subset(loo_largepool_mismatch, conc != "30 fM" & text == "*")
loo_largepool_mismatch$multiplex_label <- sub("Guides", "crRNAs", loo_largepool_mismatch$multiplex_label)

# format for plotting -----------------------------------------------------

loo_largepool_rates$conc <- factor(loo_largepool_rates$conc, levels=loo_largepool_conc)
loo_largepool_rates$mismatch <- ifelse(loo_largepool_rates$multiplex %in% 
                                         c("31 Guides", "7 Guides"),
                                       "mismatch", "no mismatch")
loo_largepool_rates$mismatch[loo_largepool_rates$multiplex == "No Protein"] <- NA
loo_largepool_rates$multiplex_label <- loo_largepool_rates$multiplex
loo_largepool_rates$multiplex_label[loo_largepool_rates$multiplex=="31 Guides"] <- "32 Guides"
loo_largepool_rates$multiplex_label[loo_largepool_rates$multiplex=="7 Guides"] <- "8 Guides"
loo_largepool_rates$multiplex_label <- factor(loo_largepool_rates$multiplex_label,
                                              levels=c("No Protein", "8 Guides", "32 Guides"))
levels(loo_largepool_rates$multiplex_label) <- c("No Protein", "8 crRNAs", "32 crRNAs")

# generate plot -----------------------------------------------------------

figure_S10 <- ggplot(subset(loo_largepool_rates, 
                                !((multiplex == "No Protein") | (conc == "30 fM"))), 
                         aes(x=factor(conc, 
                                      levels=c("RNP-only", "300 aM", "1 fM", "3 fM", "10 fM")), 
                             y=Estimate)) + 
  geom_col(aes(fill=mismatch), position=position_dodge(width=1)) + theme_classic(base_size=8) + 
  geom_errorbar(aes(ymin = Estimate + qnorm(0.025)*Std..Error,
                    ymax = Estimate - qnorm(0.025)*Std..Error,
                    fill=mismatch),
                position=position_dodge(width=1), width=0.75) + 
  facet_grid(~factor(multiplex_label, levels=c("8 crRNAs", "32 crRNAs"))) + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
  xlab("") + ylab("RFU/min") + labs(fill="") +
  scale_fill_manual(values=c("red", "grey35"),
                    labels=c("leave one\ncrRNA out", "all crRNAs")) +
  geom_text(data=loo_largepool_mismatch, aes(x=conc, y=text_height+9, label=text), size=7) + 
  geom_errorbar(data=loo_largepool_mismatch, 
                aes(ymin=text_height+7, ymax=text_height+7), size=1)

ggsave(file=file.path(figure_dir, "suppl_figure_10.pdf"),
       plot=figure_S10,
       device="pdf", width=6.5, height=3.5, units="in")
