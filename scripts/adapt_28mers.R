rm(list=ls())

library(here)
library(Biostrings)
library(ggplot2)

# load bed file of ordered guides
ordered <- read.table(file.path(here(), "outputs", "ordered_guides.bed"), 
                      col.names=c("chr", "start", "end", "id"))
ordered <- ordered[order(ordered$start, decreasing=F),]

# load results from primary screen
primary <- read.csv(file.path(here(), "outputs", "NCR_Guide_Data - Data.csv"))

# compute coordinates for ADAPT 28mers
adapt <- subset(ordered, ordered$id %in% primary$NCR.ID)
adapt$start <- adapt$start - 8

# load genome sequence
genome <- readLines(file.path(here(), "ref_data", "NC_045512v2.fa"))
genome <- genome[!grepl(">", genome)]
genome <- paste(genome, collapse="")

# pull genome sequence for ADAPT 28mers
adapt$target <- mapply(substr, start=adapt$start, stop=adapt$end, 
                       MoreArgs=list(x=genome))
write.table(adapt, quote=F, row.names=F, col.names=F, sep="\t",
            file=file.path(here(), "outputs", "adapt_28mers.tsv"))

# pull ADAPT scores
adapt_dat <- read.table(file.path(here(), "outputs", "28meroutput.tsv"),
                        header=T)
adapt_dat$target.sequence.positions <- as.numeric(sub("\\{", "", 
                                                      sub("\\}", "", 
                                                          adapt_dat$target.sequence.positions)))
adapt_dat$NCR_id <- adapt$id[match(adapt_dat$window.start+1, adapt$start)]
indices <- match(adapt_dat$NCR_id, primary$NCR.ID)
adapt_dat$activator_0fM <- primary$No.Activator.Rate[indices]
adapt_dat$activator_1pM <- primary$X1.pM.Activator.Rate[indices]
adapt_dat$activator_100fM <- primary$X100.fM.Activator.Rate[indices]
adapt_dat$activator_10fM <- primary$X10.fM.Activator.Rate[indices]
adapt_dat$activator_1fM <- primary$X1.fM.Activator.Rate[indices]
activator_concs <- c("0fM", "1fM", "10fM", "100fM", "1pM")

# plot ADAPT median activity by rate
median_activity <- lapply(activator_concs,
                          function(x) {
                            data.frame(ADAPT=adapt_dat$guide.set.median.activity,
                                       screen=as.numeric(adapt_dat[,paste0("activator_", x)]),
                                       activator_conc=x)
                          })
median_activity <- do.call(rbind, median_activity)
median_activity$activator_conc <- factor(median_activity$activator_conc, levels=activator_concs)
ggplot(median_activity, aes(x=ADAPT, y=screen, col=activator_conc)) +
  geom_point() + geom_smooth(method="lm", formula=y~x) +
  facet_grid(activator_conc ~ ., scales="free_y") + 
  theme_bw() + xlab("ADAPT: guide set median activity") + ylab("rate from primary screen") + 
  labs(col="activator concentration")

# plot ADAPT rank by rate rank
median_rank <- lapply(activator_concs,
                        function(x) {
                          data.frame(ADAPT=order(adapt_dat$guide.set.median.activity,
                                                 decreasing=T),
                                     screen=order(as.numeric(adapt_dat[,paste0("activator_", x)]),
                                                  decreasing=T),
                                     activator_conc=x)
                        })
median_rank <- do.call(rbind, median_rank)
median_rank$activator_conc <- factor(median_activity$activator_conc, levels=activator_concs)
ggplot(median_rank, aes(x=ADAPT, y=screen, col=activator_conc)) +
  geom_point() + geom_smooth(method="lm", formula=y~x) +
  facet_grid(activator_conc ~ .) + 
  theme_bw() + xlab("ADAPT: guide set median activity (rank)") + ylab("rate from primary screen (rank)") + 
  labs(col="activator concentration")

# combine plots
plot_dat <- rbind(data.frame(median_activity, type="activity"),
                  data.frame(median_rank, type="rank"))
ggplot(plot_dat, aes(x=ADAPT, y=screen, col=activator_conc)) +
  geom_point() + geom_smooth(method="lm", formula=y~x) +
  facet_wrap(facets=c("activator_conc", "type"), scales="free", ncol=2) +
  theme_bw() + xlab("ADAPT: median") + ylab("rate from primary screen") +
  labs(col="activator concentration")

# get 28mer guide sequences
adapt12_fname <- file.path(here(), "outputs", "AdaptSelections.xlsx")
adapt12 <- xlsx::read.xlsx(adapt12_fname, sheetIndex=1)
adapt12$guide <- paste0("uagaccaccccaaaaaugaaggggacuaaaac",
                        as.character(RNAStringSet(reverseComplement(DNAStringSet(adapt12$target.sequences)))))
xlsx::write.xlsx(adapt12, adapt12_fname)
