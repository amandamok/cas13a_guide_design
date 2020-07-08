rm(list=ls())

library(here())

# load data
plusStrand <- read.table(file.path(here(), "outputs", "cas13a_20nt", "cas13a_results_summary.txt"))
minusStrand <- read.table(file.path(here(), "outputs", "cas13a_minusStrand_20nt", "cas13a_minusStrand_results_summary.txt"))
dat <- rbind(plusStrand, minusStrand)
offtarget_vir <- read.csv(file.path(here(), "ref_data", "offtarget_analysis", "20200615_1056_guides.csv"))
offtarget_bac <- read.table(file.path(here(), "ref_data", "offtarget_analysis", "mat_refbac_20200609.tsv"), 
                                      header=T, stringsAsFactors=F)
offtarget_bac$hits_01 <- offtarget_bac$X0 + offtarget_bac$X1

# add offtarget data
dat$offtarget_vir <- NA
# dat$offtarget_vir[prodlim::row.match(offtarget_candidates[, c("start", "strand")], dat[, c("start", "strand")])] <- T
dat$offtarget_vir[prodlim::row.match(offtarget_vir[, c("start", "strand")], dat[, c("start", "strand")])] <- F
dat$offtarget_bac <- NA
dat$offtarget_bac[prodlim::row.match(offtarget_bac[, c("start_position", "strand")], 
                                     dat[, c("start", "strand")])] <- (offtarget_bac$hits_01>0)
dat$offtarget_any <- !(!dat$offtarget_bac & !dat$offtarget_vir)
offtarget_dat <- subset(dat, !offtarget_any)

# load amplicons
amplicons <- readLines(file.path(here(), "isothermal_amplification", "winner_amplicons_mapped.sam"))
amplicons <- amplicons[!grepl("^@", amplicons)]
amplicons <- data.frame(matrix(unlist(strsplit(amplicons, split="\t")), nrow=length(amplicons), byrow=T), stringsAsFactors=F)
amplicons <- amplicons[,c(1, 4, 10)]
colnames(amplicons) <- c("amplicon", "pos", "seq")
amplicons$pos <- as.numeric(as.character(amplicons$pos))
amplicons$end <- amplicons$pos + nchar(amplicons$seq)

# write bed files

fwd_strand <- subset(offtarget_dat, strand=="+")
fwd_strand <- data.frame(chrom = "NC_045512v2",
                         chromStart = fwd_strand$start,
                         chromEnd = fwd_strand$start+20-1,
                         name = paste("+", fwd_strand$start, sep="_"))
write.table(fwd_strand, file=file.path(here(), "outputs", "guides_fwd.bed"), quote=F, sep="\t",
            row.names=F, col.names=F)

rev_strand <- subset(offtarget_dat, strand=="-")
rev_strand <- data.frame(chrom = "NC_045512v2",
                         chromStart = rev_strand$start,
                         chromEnd = rev_strand$start+20-1,
                         name = paste("+", rev_strand$start, sep="_"))
write.table(rev_strand, file=file.path(here(), "outputs", "guides_rev.bed"), quote=F, sep="\t",
            row.names=F, col.names=F)

amplicons <- data.frame(chrom = "NC_045512v2",
                        chromStart = amplicons$pos,
                        chromEnd = amplicons$end,
                        name = amplicons$amplicon)
write.table(amplicons, file=file.path(here(), "outputs", "amplicons.bed"), quote=F, sep="\t",
            row.names=F, col.names=F)

# figure out non-overlapping regions
fwd_strand$region <- NA
for(x in seq(nrow(fwd_strand))) {
  if(x == 1) {
    region <- 1
    region_start <- fwd_strand$start[x]
  } else {
    if(fwd_strand$start[x] > region_start + 20 + 10) {
      region <- region + 1
      region_start <- fwd_strand$start[x]
    }
  }
  fwd_strand$region[x] <- region
}
rev_strand$region <- NA
for(x in seq(nrow(rev_strand))) {
  if(x == 1) {
    region <- 1
    region_start <- rev_strand$chromStart[x]
  } else {
    if(rev_strand$chromStart[x] > region_start + 20 + 10) {
      region <- region + 1
      region_start <- rev_strand$chromStart[x]
    }
  }
  rev_strand$region[x] <- region
}
