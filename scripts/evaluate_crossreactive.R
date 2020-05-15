##################################################
### EUA cross-reactivity

# in silico cross-reactivity (FDA): >80% sequence homology between one of the primers/probes
# and any sequence present in the targeted microorganism

# recommended list:
## other high priority pathogens from the same genetic family
# - human coronavirus 229E: NC_002645.1
# - human coronavirus OC43: NC_006213.1
# - human coronavirus HKU1: NC_006577.2
# - human coronavirus NL63: NC_005831.2
# - SARS-coronavirus: NC_004718.3
# - MERS-coronavirus: NC_019843.3
## high priority organisms likely in circulating areas
# - adenovirus (e.g. C1 Ad. 71)
# - human metapneumovirus (hMPV)
# - parainfluenza virus 1-4
# - influenza A & B
# - enterovirus (e.g. EV68)
# - respiratory syncytial virus
# - rhinovirus
# - chlamydia pneumoniae
# - haemophilus influenzae
# - legionella pneumophila
# - mycobacterium tuberculosis
# - streptococcus pneumoniae
# - streptococcus pyogenes
# - bordetella pertussis
# - mycoplasma pneumoniae
# - pneumocystis jirovecii (PJP)
# - candida albicans
# - pseudomonas aeruginosa
# - staphylococcus epidermidis
# - staphylococcus salivarius
# - pooled human nasal wash

library(optparse)
library(Biostrings)
library(ggplot2)

option_list <- list(make_option(c("-m", "--mismatch"), type="integer", default=4,
                                help="number of mismatches allowed", metavar="integer"),
                    make_option(c("-t", "--target"), type="character", default=NULL,
                                help="target is DNA or RNA", metavar="DNA or RNA"),
                    make_option(c("-o", "--out"), type="character", default=".", 
                                help="output directory", metavar="character")) 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(!("target") %in% names(opt)) {
  cat("\nERROR: no target specified\n")
  q(save="no")
} else {
  cat("\nEvaluating cross-reactivity\n")
}

# read in probes
windows <- read.table(file.path(opt$out, "windows.txt"), header=T, stringsAsFactors=F)
targets <- suppressWarnings(PDict(windows$target, max.mismatch=opt$mismatch))

# read in off-target genomes
human_CoV_dir <- "~/covid-19/ref_data/human_CoV/"
human_CoVs <- c("huCoV_229E", "SARS", "huCoV_NL63", "huCoV_OC43", "huCoV_HKU1", "MERS")
other_pathogens_dir <- "~/covid-19/ref_data/other_pathogens/"
offtarget_ssRNA_plusStrand <- c("enterovirus", "rhinovirus", "parechovirus")
offtarget_ssRNA_minusStrand <- c("hMPV", paste0("parainfluenza_", 1:4), paste0("influenza_", c("A", "B")), "resp_syncytial_A")
offtarget_dsDNA <- c("adenovirus", "c_pneumoniae", "h_influenzae", "l_pneumophila", "m_tuberculosis", "s_pneumoniae", "s_pyogenes",
                     "b_pertussis", "m_pneumoniae", "p_jirovecii", "c_albicans", "p_aeruginosa", "s_epidermidis", "s_salivarius")
if(opt$target == "DNA") {
  offtargets_fname <- c(paste0(human_CoV_dir, list.files(human_CoV_dir, pattern="*fa")),
                        paste0(other_pathogens_dir, c(offtarget_ssRNA_plusStrand, offtarget_ssRNA_minusStrand, offtarget_dsDNA), ".fa"))
} else {
  offtargets_fname <- c(paste0(human_CoV_dir, list.files(human_CoV_dir, pattern="*fa")),
                        paste0(other_pathogens_dir, offtarget_ssRNA_plusStrand, ".fa"),
                        paste0(other_pathogens_dir, offtarget_ssRNA_minusStrand, ".fa"),
                        paste0(other_pathogens_dir, offtarget_dsDNA, "_rna.fa"))
}
offtargets <- lapply(offtargets_fname, readDNAStringSet)
names(offtargets) <- c(human_CoVs, offtarget_ssRNA_plusStrand, offtarget_ssRNA_minusStrand, offtarget_dsDNA)

# count number of alignments, allowing for 4 mismatches
alignment_cts <- sapply(seq_along(offtargets),
                        function(x) {
                          tmp_offtarget <- offtargets[[x]]
                          cts <- rowSums(sapply(seq_along(tmp_offtarget),
                                                function(y) {
                                                  countPDict(targets, unlist(tmp_offtarget[y]), max.mismatch=opt$mismatch)
                                                }))
                          if((opt$target=="DNA") & (names(offtargets)[x] %in% offtarget_dsDNA)) {
                            tmp_offtarget <- reverseComplement(tmp_offtarget)
                            cts <- cts + rowSums(sapply(seq_along(tmp_offtarget),
                                                        function(y) {
                                                          countPDict(targets, unlist(tmp_offtarget[y]), max.mismatch=opt$mismatch)
                                                        }))
                          }
                          return(cts)
                        })
colnames(alignment_cts) <- names(offtargets)

# compile results table
results <- data.frame(start=windows$start,
                      strand=windows$strand,
                      offtarget_rate=rowMeans(alignment_cts>0),
                      offtargets=sapply(seq.int(nrow(alignment_cts)),
                                        function(x) {
                                          paste(colnames(alignment_cts)[alignment_cts[x,]>0], collapse=",")
                                        }))

# output off-target alignments
write.table(data.frame(start=windows$start,
                       strand=windows$strand,
                       alignment_cts),
            file=file.path(opt$out, paste0("crossreactive_offtargets_", opt$target, ".txt")),
            quote=F, sep="\t", row.names=F)
write.table(results,
            file=file.path(opt$out, paste0("score_crossreactive_", opt$target, ".txt")),
            quote=F, sep="\t", row.names=F)

# plot results
alignment_summary <- data.frame(target=factor(colnames(alignment_cts), levels=colnames(alignment_cts)),
                                offtarget_rate=colMeans(alignment_cts>0),
                                target_type=NA)
alignment_summary$target_type[colnames(alignment_cts) %in% human_CoVs] <- "human CoVs"
alignment_summary$target_type[colnames(alignment_cts) %in% offtarget_ssRNA_minusStrand] <- "ssRNA (-)"
alignment_summary$target_type[colnames(alignment_cts) %in% offtarget_ssRNA_plusStrand] <- "ssRNA (+)"
alignment_summary$target_type[colnames(alignment_cts) %in% offtarget_dsDNA] <- "dsDNA"
crossreactive_plot <- ggplot(alignment_summary, aes(target, offtarget_rate)) + 
  geom_col(aes(fill=target_type)) + theme_bw() + theme(axis.text.x=element_text(angle=90)) + 
  ggtitle(paste0("off-target cross-reactivity: ", opt$target), 
          subtitle=paste0("all spacers (n=", nrow(alignment_cts), ")")) + 
  xlab("") + ylab("% spacers aligned with ≤4 mismatches")
ggsave(filename=file.path(opt$out, paste0("crossreactive_offtargets_", opt$target, ".txt")),
       crossreactive_plot, device="pdf", width=6, height=4, units="in")

# ordered <- readLines("~/covid-19/outputs/spacers_20200504.txt")
# ordered <- gsub("T", "U", ordered)
# alignment_ordered <- alignment_cts[match(ordered, windows$spacer),]
# alignment_ordered_summary <- alignment_summary
# alignment_ordered_summary$offtarget_rate <- colMeans(alignment_ordered>0)
# ggplot(alignment_ordered_summary, aes(target, offtarget_rate)) +
#   geom_col(aes(fill=target_type)) + theme_bw() + theme(axis.text.x=element_text(angle=90)) +
#   ggtitle(paste0("off-target cross-reactivity: ", opt$target), 
#           subtitle=paste0("all spacers (n=", nrow(alignment_ordered_summary), ")")) + 
#   xlab("") + ylab("% spacers aligned with ≤4 mismatches")
