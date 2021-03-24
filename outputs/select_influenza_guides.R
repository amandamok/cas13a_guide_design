rm(list=ls())

library(here)
library(ggplot2)

output_dir <- file.path(here(), "outputs")

# influenza A (H1N1)
A_H1N1 <- read.table(file.path(output_dir, "influenzaA_H1N1_cas13a_neg_20nt",
                               "influenzaA_H1N1_cas13a_neg_results_summary.txt"))
A_H1N1$virus <- "influenzaA_H1N1"
# influenza A (H3N2)
A_H3N2 <- read.table(file.path(output_dir, "influenzaA_H3N2_cas13a_neg_20nt",
                               "influenzaA_H3N2_cas13a_neg_results_summary.txt"))
A_H3N2$virus <- "influenzaA_H3N2"
# influenza B
B <- read.table(file.path(output_dir, "influenzaB_human_2015-2020_cas13a_neg_20nt",
                          "influenzaB_cas13a_neg_results_summary.txt"))
B$virus <- "influenzaB"

# combine guides (n=40694)
all_guides <- rbind(A_H1N1, A_H3N2, B)
# restrict to good secondary structure (n=8212)
selected_guides <- subset(all_guides, has_crRNA_hairpin & crRNA_spacer_basepairs <=4)
# restrict to sensitivity ≥ 90% (n=6654)
selected_guides <- subset(selected_guides, sensitivity>=0.95)
# annotate all_guides
all_guides$selected <- "filtered"
all_guides$selected[which(with(all_guides,
                               has_crRNA_hairpin &
                                 crRNA_spacer_basepairs <= 4 &
                                 sensitivity >= 0.95))] <- "selected"
all_guides$selected <- factor(all_guides$selected, levels=c("selected", "filtered"))

# plots
all_guides_sensitivity <- ggplot(all_guides, aes(x=start, y=sensitivity, col=selected)) +
  geom_point(alpha=0.5, size=0.5) + facet_grid(virus ~ segment, scales="free_x") +
  ggtitle("all guides") + theme_bw() + scale_color_manual(values=c("red", "black")) +
  xlab("genomic position") + ylab("sensitivity (allowing 1 mismatch)")
all_guides_structure <- ggplot(all_guides, aes(x=has_crRNA_hairpin,
                                               y=crRNA_spacer_basepairs,
                                               col=selected)) +
  geom_jitter(alpha=0.5, size=0.5) + facet_grid(virus ~ segment) + theme_bw() +
  geom_label(data=subset(aggregate(crRNA_spacer_basepairs ~ virus + segment + has_crRNA_hairpin,
                                   data=all_guides, FUN=length),
                         has_crRNA_hairpin),
             aes(x=has_crRNA_hairpin, y=0, label=crRNA_spacer_basepairs), col="red") +
  xlab("intact crRNA hairpin structure") + ylab("number basepairs in spacer") +
  ggtitle("all guides")
ggplot(selected_guides, aes(x=start, y=sensitivity)) + geom_point(alpha=0.2, size=0.75) +
  facet_grid(virus ~ segment, scales="free_x") + theme_bw() +
  xlab("genomic position") + ylab("sensitivity (allowing 1 mismatch)")

# generic guides
duplicated_guides <- selected_guides$spacer[duplicated(selected_guides$spacer)]
duplicated_guides <- subset(selected_guides, spacer %in% duplicated_guides)
selected_guides$target_virus <- selected_guides$virus
selected_guides$target_virus[selected_guides$spacer %in% duplicated_guides$spacer] <- "influenzaA_H1N1_H3N2"

# write output
selected_guides$guide_sequence <- paste0("UAGACCACCCCAAAAAUGAAGGGGACUAAAAC",
                                         selected_guides$spacer)
write.table(selected_guides, file=file.path(output_dir, "influenza_guides_20210324.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T)
