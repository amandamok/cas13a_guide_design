rm(list=ls())

library(prodlim)
library(here)

setwd(here())

# load all guides
plusStrand <- read.table("outputs/cas13a_20nt/cas13a_results_summary.txt")
minusStrand <- read.table("outputs/cas13a_minusStrand_20nt/cas13a_minusStrand_results_summary.txt")
dat <- rbind(plusStrand, minusStrand)

# load offtarget data
offtarget_candidates <- read.table("outputs/candidate_spacers.txt", stringsAsFactors=F, header=T)
offtarget_vir <- read.csv("ref_data/offtarget_analysis/20200615_1056_guides.csv")
offtarget_bac <- read.table("ref_data/offtarget_analysis/mat_refbac_20200609.tsv", header=T, stringsAsFactors=F)
offtarget_bac$hits_01 <- offtarget_bac$X0 + offtarget_bac$X1

# load annotation regions
shape_map <- read.table("ref_data/shape_map_ssRNA.bed", col.names=c("chr", "start", "stop"))
shape_map <- subset(shape_map, stop-start+1 >= 20)
dms_map <- read.table("ref_data/dms_map_unstructured.bed", col.names=c("chr", "start", "stop"))
genome <- data.frame(pos=1:29903,
                     shape=F,
                     dms=F)
for(x in seq(nrow(shape_map))) { genome$shape[shape_map$start[x]:shape_map$stop[x]] <- T }
for(x in seq(nrow(dms_map))) { genome$dms[dms_map$start[x]:dms_map$stop[x]] <- T}
dms_map <- subset(dms_map, stop-start+1 >= 20)
gblocks <- read.table("ref_data/gblocks.bed", col.names=c("chr", "start", "stop", "gblock"))

# structure annotations
shape_map_overlap <- do.call(rbind,
                             lapply(seq(nrow(shape_map)),
                                    function(x) {
                                      tmp <- subset(dat, start>=shape_map$start[x] & start<=(shape_map$stop[x]-20+1))
                                      if(nrow(tmp)==0) { return(NULL) }
                                      else { return(tmp) }
                                    }))
dat$shape_map <- F
dat$shape_map[row.match(shape_map_overlap[, c("start", "strand")], dat[, c("start", "strand")])] <- T
dms_map_overlap <- do.call(rbind,
                           lapply(seq(nrow(dms_map)),
                                  function(x) {
                                    tmp <- subset(dat, start>=dms_map$start[x] & start<=(dms_map$stop[x]-20+1))
                                    if(nrow(tmp)==0) { return(NULL) }
                                    else { return(tmp) }
                                  }))
dat$dms_map <- F
dat$dms_map[row.match(dms_map_overlap[, c("start", "strand")], dat[, c("start", "strand")])] <- T

# gblock annotations
dat$gblock <- NA
gblock_margin <- 30
for(x in seq(nrow(gblocks))) {
  which_guides <- (dat$start > (gblocks$start[x]+gblock_margin)) & ((dat$start+20) < (gblocks$stop[x]-gblock_margin))
  dat$gblock[which_guides] <- gblocks$gblock[x]
}

# subset to candidates, annotate with offtarget hits
candidates <- dat[row.match(offtarget_candidates[, c("start", "strand")], dat[, c("start", "strand")]),]
candidates$no_offtarget_vir <- F
candidates$no_offtarget_vir[row.match(offtarget_vir[, c("start", "strand")],
                                      candidates[, c("start", "strand")])] <- T
candidates$no_offtarget_bac <- (offtarget_bac$hits_01[row.match(candidates[, c("start", "strand")],
                                                                offtarget_bac[, c("start_position", "strand")])]==0)


# candidates: shape_map & dms_map
candidates_yesMap_yesDms <- subset(candidates, shape_map & dms_map & strand=="+")
table(candidates_yesMap_yesDms$no_offtarget_bac, candidates_yesMap_yesDms$no_offtarget_vir)

# candidates: shape_map & !dms_map
candidates_yesMap_noDms <- subset(candidates, shape_map & !dms_map & strand=="+")
table(candidates_yesMap_noDms$no_offtarget_bac, candidates_yesMap_noDms$no_offtarget_vir)

# candidates: !shape_map & dms_map
candidates_noMap_yesDms <- subset(candidates, !shape_map & dms_map & strand=="+")
table(candidates_noMap_yesDms$no_offtarget_bac, candidates_noMap_yesDms$no_offtarget_vir)

# candidates: no offtarget hits for bacteria and viruses
candidates_noOfftarget <- subset(candidates, no_offtarget_bac & no_offtarget_vir & strand=="+")

# select 96 to order
selection <- unique(rbind(candidates_yesMap_yesDms,
                          subset(candidates_yesMap_noDms, no_offtarget_vir | no_offtarget_bac),
                          subset(candidates_noMap_yesDms, no_offtarget_vir | no_offtarget_bac),
                          candidates_noOfftarget))
selection <- subset(selection, !is.na(gblock))
selection <- selection[order(selection$start), ]
selection$region <- NA
for(x in seq(nrow(selection))) {
  if(x == 1) {
    region <- 1
    region_start <- selection$start[x]
  } else {
    if(selection$start[x] > region_start+30) {
      region <- region + 1
      region_start <- selection$start[x]
    }
  }
  selection$region[x] <- region
}

write.csv(selection, file="outputs/guide_selection_20200708.csv", quote=F, row.names=F)

# ott guides
ott_guides <- readLines("ref_data/guides_ott.fa")
ott_guides <- data.frame(id = sub("> ", "", ott_guides[grepl(">", ott_guides)]),
                         spacer = ott_guides[!grepl(">", ott_guides)], 
                         stringsAsFactors=F)
ott_guides_subset <- subset(dat, spacer %in% ott_guides$spacer & strand == "+")
ott_guides_subset$id <- ott_guides$id[match(ott_guides$spacer, ott_guides_subset$spacer)]
write.csv(ott_guides_subset, file="outputs/ott_guides.csv", quote=F, row.names=F)
