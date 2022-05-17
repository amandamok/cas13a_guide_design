rm(list=ls())

library(here)
library(ggplot2)

gblocks_bed <- read.table(file.path(here(), "ref_data", "gblocks.bed"),
                          col.names=c("seq", "start", "end", "gblock"))
gblocks_rate <- read.csv(file.path(here(), "NCR", "data", "gblock_rate.csv"))
guide_features <- read.csv(file.path(here(), "NCR", "data", "guide_features.csv"))

# add guide coordinate to gblock rates
gblocks_rate$start <- guide_features$start[match(gblocks_rate$guide_id,
                                                 sub("NCR_", "", guide_features$NCR.id))]

# add antitag label to gblock rates
gblocks_rate$antitag <- guide_features$antitag_pos1[match(gblocks_rate$guide_id,
                                                          sub("NCR_", "", guide_features$NCR.id))]

# annotate guides by gblock
gblocks_rate$gblock <- NA
for(x in seq(nrow(gblocks_bed))) {
  gblocks_rate$gblock[(gblocks_rate$start > gblocks_bed$start[x]) & 
                        ((gblocks_rate$start + 19) < gblocks_bed$end[x])] <- gblocks_bed$gblock[x]
}

# remove outlier: NCR_13552

# plot guide rates by gblock
ggplot(subset(gblocks_rate, guide_id != 1352), 
       aes(x=paste0("NCR_", guide_id), y=Estimate, fill=antitag,
           ymin=Estimate + qnorm(0.025)*Std.Error,
           ymax=Estimate + qnorm(0.975)*Std.Error,)) + 
  geom_col() + geom_errorbar() + theme_classic() + 
  xlab("") + ylab("rate above RNP-only [RFU/min]") + 
  ggtitle("", subtitle="activator: 100 fM pooled gBlocks") + 
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5)) + 
  facet_wrap(~gblock, scales="free_x")
