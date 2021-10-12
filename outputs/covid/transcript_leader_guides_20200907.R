#####
# design 10 guides that target 5' transcript leader (pos 1-100)

rm(list=ls())

library(here)

plusStrand <- read.table(file.path(here(), "outputs", "cas13a_20nt", 
                                   "cas13a_results_summary.txt"),
                         stringsAsFactors=F)
minusStrand <- read.table(file.path(here(), "outputs", "cas13a_minusStrand_20nt", 
                                    "cas13a_minusStrand_results_summary.txt"),
                          stringsAsFactors=F)
dat <- rbind(plusStrand, minusStrand, stringsAsFactors=F)

# 1. subset to target region
dat <- subset(dat, start <= 100 - unique(nchar(dat$spacer)))
# 2. subset to high sensitivity
dat <- subset(dat, sensitivity_01 > 0.9)
# 3. subset to good secondary structure
dat <- subset(dat, has_crRNA_hairpin)
dat <- subset(dat, crRNA_spacer_basepairs <= 4)
# 4. subset to forward strand
dat <- subset(dat, strand=="+")

# add complete crRNA sequence
crRNA_repeat <- toupper("uagaccaccccaaaaaugaaggggacuaaaac")
dat$crRNA_sequence <- paste0(crRNA_repeat, dat$spacer)

# select 10 guides: remove guide w/ highest cross-reactivity against human transcriptome
dat <- dat[-which.max(dat$match_against_hg38), ]

write.table(dat, 
            file=file.path(here(), "outputs", "transcript_leader_guides_20200907.txt"),
            quote=F, sep="\t", row.names=F, col.names=T)
