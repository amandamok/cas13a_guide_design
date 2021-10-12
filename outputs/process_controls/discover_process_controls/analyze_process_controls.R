rm(list=ls())

rnafold <- readLines("crRNAs_RNAfold.txt")
num_crRNAs <- length(rnafold)/2
direct_repeat <- "UAGACCACCCCAAAAAUGAAGGGGACUAAAAC"

rnafold <- data.frame(guide=paste0(rep(c("ACTB", "GAPDH9", "GAPDH5"), each=3),
                                   rep(seq(3), times=3)),
                      crRNA = rnafold[2*seq(num_crRNAs)-1],
                      crRNA_structure = sapply(rnafold[2*seq(num_crRNAs)],
                                               function(x) {
                                                 strsplit(x, split=" ")[[1]][1]
                                               }),
                      row.names=NULL)
rnafold$spacer <- sapply(as.character(rnafold$crRNA),
                         function(x) {
                           sub(direct_repeat, "", x)
                         })
rnafold$dr_structure <- sapply(as.character(rnafold$crRNA_structure),
                               function(x) {
                                 substr(x, 1, nchar(direct_repeat))
                               })
rnafold$spacer_structure <- sapply(as.character(rnafold$crRNA_structure),
                                   function(x) {
                                     substr(x, nchar(direct_repeat)+1, nchar(x))
                                   })
rnafold$num_bp_spacer <- sapply(rnafold$spacer_structure,
                                function(x) {
                                  tmp_structure <- strsplit(x, split="")[[1]]
                                  return(sum(tmp_structure %in% c("(", ")")))
                                })

write.csv(rnafold, file="process_control_results.csv", row.names=F)
