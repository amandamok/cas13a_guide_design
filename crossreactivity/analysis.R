rm(list=ls())

library(here)

parse_ct <- function(ct_fname) {
  # report minimum folding energy or maximum number basepairs
  ## ct_fname: character ; path to .ct output file from DuplexFold (RNAstructure)
  dat <- readLines(ct_fname)
  which_energy <- grep("ENERGY", dat)
  if(length(which_energy)==0) {
    min_energy <- NA
    which_energy <- 1
  } else {
    dat_energy <- sapply(which_energy,
                         function(x) {
                           as.numeric(strsplit(dat[x], split=" ")[[1]][8])
                         })
    min_energy <- min(dat_energy)
  }
  dat_bp <- sapply(which_energy,
                   function(x) {
                     tmp_bp <- dat[(x+1):(x+43)]
                     tmp_bp <- sapply(tmp_bp,
                                      function(x) {
                                        tmp <- strsplit(x, split=" ")[[1]]
                                        tmp <- tmp[tmp!=""]
                                        return(as.numeric(tmp[5]))
                                      })
                     return(sum(tmp_bp!=0)/2)
                   })
  # report outputs
  return(list(min_energy=min_energy, max_bp=max(dat_bp)))
}

# write guide sequences to individual fasta files
guides <- read.csv(file.path(here(), "outputs", "NCR_Guide_Data - Data.csv"),
                   stringsAsFactors=F)[, c("NCR.ID", "Sequence")]
guide_fasta_dir <- file.path(here(), "crossreactivity", "guides")
if(!dir.exists(file.path(guide_fasta_dir))) { dir.create(file.path(guide_fasta_dir)) }
for(x in seq(nrow(guides))) {
  tmp_ID <- guides$NCR.ID[x]
  tmp_seq <- sub("uagaccaccccaaaaaugaaggggacuaaaac", "", guides$Sequence[x])
  tmp_guide <-
    writeLines(c(paste0(">", tmp_ID), tmp_seq),
               con=file.path(guide_fasta_dir, paste0(tmp_ID, ".fa")))
}

# write script to calculate paired structures
structures_dir <- file.path(here(), "crossreactivity", "structures")
if(!dir.exists(structures_dir)) { dir.create(structures_dir) }
dat <- data.frame(t(combn(guides$NCR.ID, m=2)))
colnames(dat) <- c("guide_1", "guide_2")
duplexfold_commands <- sapply(seq(nrow(dat)),
                              function(x) {
                                paste("DuplexFold",
                                      file.path(guide_fasta_dir, paste0(dat$guide_1[x], ".fa")),
                                      file.path(guide_fasta_dir, paste0(dat$guide_2[x], ".fa")),
                                      file.path(structures_dir, paste0(dat$guide_1[x], "_", dat$guide_2[x], ".ct")))
                              })
duplexfold_commands <- c("#!/bin/bash", "", duplexfold_commands)
writeLines(duplexfold_commands, con=file.path("crossreactivity", "duplexfold_commands.sh"))
system(paste("chmod +x", file.path("crossreactivity", "duplexfold_commands.sh")))
# run duplexfold_commands.sh in separate terminal outside RStudio

# parse folding output
dat$min_energy <- NA
dat$max_bp <- NA
for(x in seq(nrow(dat))) {
  tmp_fname <- file.path(structures_dir, paste0(dat$guide_1[x], "_", dat$guide_2[x], ".ct"))
  if(file.exists(tmp_fname)) {
    tmp_output <- parse_ct(tmp_fname)
    dat$min_energy[x] <- tmp_output$min_energy
    dat$max_bp[x] <- tmp_output$max_bp
  }
}

# plot minimum free energy
plot_min_energy <- ggplot(dat, aes(x=guide_1, y=guide_2, fill=min_energy)) + geom_tile(col=1) + theme_bw() +
  scale_fill_gradient2(low="navyblue", mid="indianred1", na.value="white") + xlab("") + ylab("") +
  ggtitle("Minimum free energy (RNAstructure::DuplexFold)",
          subtitle=paste("min =", min(dat$min_energy, na.rm=T), ", max =", max(dat$min_energy, na.rm=T))) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), panel.grid.major.y=element_blank())
ggsave(plot_min_energy, filename=file.path("crossreactivity", "plot_min_energy.pdf"),
       device="pdf", width=12, height=12, units="in", dpi="print")

# plot maximum number base-pairing
plot_max_bp <- ggplot(dat, aes(x=guide_1, y=guide_2, fill=max_bp)) + geom_tile(col=1) + theme_bw() +
  scale_fill_gradient2(low="navyblue", mid="greenyellow", na.value="white") + xlab("") + ylab("") +
  ggtitle("Maximum number of base-paired positions (RNAstructure::DuplexFold)",
          subtitle=paste("min =", min(dat$max_bp, na.rm=T), ", max =", max(dat$max_bp, na.rm=T))) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), panel.grid.major.y=element_blank())
ggsave(plot_max_bp, filename=file.path("crossreactivity", "plot_min_energy.pdf"),
       device="pdf", width=12, height=12, units="in", dpi="print")
