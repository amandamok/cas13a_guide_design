rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(ggseqlogo)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# functions ---------------------------------------------------------------

compute_nt_freq <- function(sequences, bases=c("A", "U", "C", "G")) {
  # sequences: character vector; sequences of same length
  # bases: character vector; letters to compute frequences of at each position
  seq_length <- unique(nchar(sequences))
  seq_nt <- matrix(unlist(strsplit(sequences, split="")), ncol=seq_length, byrow=T)
  seq_freq <- sapply(seq(seq_length),
                     function(x) {
                       sapply(bases,
                              function(base) {
                                mean(seq_nt[,x]==base)
                              })
                     })
  rownames(seq_freq) <- bases
  return(data.frame(t(seq_freq)))
}

calculate_relative_info <- function(q_ik, p_ik) {
  # q_ik: data.frame; base frequencies per position in queries (bases as columns)
  # p_ik: data.frame; base frequencies per position in background (bases as columns)
  num_positions <- nrow(q_ik)
  if(nrow(p_ik) != num_positions) {
    stop("q_ik and p_ik have differing number of positions")
  }
  if(any(!(colnames(q_ik) %in% colnames(p_ik)))) {
    stop("q_ik and p_ik have differing alphabets")
  }
  info_content <- sapply(seq(num_positions),
                         function(x) {
                           sum(q_ik[x,] * log2(q_ik[x,] / p_ik[x,]))
                         })
  denominator <- sapply(seq(num_positions),
                        function(x) {
                          sum(q_ik[x,] / p_ik[x,])
                        })
  base_heights <- sapply(colnames(q_ik),
                         function(base) {
                           q_ik[[base]]/p_ik[[base]]/denominator*info_content
                         })
  return(t(base_heights))
}

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))

# calculate nt frequencies ------------------------------------------------

num_spacers <- 50

# nt frequencies in all spacers
allSpacers_freq <- compute_nt_freq(subset(guide_rate$spacer,
                                          nchar(guide_rate$spacer)==20))
nt_freq <- sapply(c("A", "U", "C", "G"),
                  function(x) mean(unlist(strsplit(subset(guide_rate$spacer,
                                                          nchar(guide_rate$spacer)==20), 
                                                   split=""))==x))

# top spacers
topSpacers_freq <- compute_nt_freq(subset(guide_rate$spacer,
                                          nchar(guide_rate$spacer) == 20 &
                                            rank(-guide_rate$Estimate) <= num_spacers))
topSpacers_info <- calculate_relative_info(topSpacers_freq, allSpacers_freq)

# bottom spacers
bottomSpacers_freq <- compute_nt_freq(subset(guide_rate$spacer,
                                             nchar(guide_rate$spacer) == 20 &
                                               rank(guide_rate$Estimate) <= num_spacers))
bottomSpacers_info <- calculate_relative_info(bottomSpacers_freq, allSpacers_freq)

# generate plot -----------------------------------------------------------

figure_S4B <- (ggseqlogo(topSpacers_info, method="custom") + 
   ggtitle(paste("Top", num_spacers, "spacers")) + ylim(0, 0.5) + 
   theme_classic(base_size=8) + ylab("bits")) / 
  (ggseqlogo(bottomSpacers_info, method="custom") +
     ggtitle(paste("Bottom", num_spacers, "spacers")) + ylim(0, 0.5) +
     theme_classic(base_size=8) + ylab("bits") + 
     xlab("position relative to spacer"))

ggsave(file=file.path(figure_dir, "suppl_figure_4B.pdf"),
       plot=figure_S4B,
       device="pdf", width=6.5, height=3, units="in")
