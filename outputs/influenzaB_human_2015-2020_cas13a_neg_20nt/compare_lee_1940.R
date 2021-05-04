library(here)
library(Biostrings)
library(ggplot2)
library(patchwork)

B_ref_dir <- file.path(here(), "ref_data", "influenza")
fasta_files <- grep("consensus", invert=T, value=T,
                    list.files(file.path(B_ref_dir, "influenzaB_human_2015-2021")))
consensus_fname <- grep("allSegments", value=T,
                        list.files(file.path(B_ref_dir, "influenzaB_human_2015-2021")))

B_Lee <- readDNAStringSet(file.path(B_ref_dir, "B_Lee_1940.fasta"))
B_Victoria <- readDNAStringSet(file.path(B_ref_dir, "B_Victoria_504_2000.fasta"))
B_consensus <- readDNAStringSet(file.path(B_ref_dir, "influenzaB_human_2015-2021",
                                          consensus_fname))

calculate_match <- function(sample_segments, reference_segment) {
  # sample_segments: DNAStringSet; segment sequences of samples
  # reference_segment: DNAStringSet; segment sequence of reference
  # process inputs
  reference_segment_string <- unlist(strsplit(as.character(reference_segment), split=""))
  sample_segments <- gsub("-", "", sample_segments)
  # align samples to reference
  alignment <- as.matrix(pairwiseAlignment(sample_segments, reference_segment))
  # calculate percent match
  percent_match <- sapply(seq(nrow(alignment)),
                                 function(x) {
                                   tmp_length <- width(sample_segments)[x]
                                   tmp_match <- sum(alignment[x,] == reference_segment_string)
                                   return(tmp_match/tmp_length)
                                 })
  # calculate positional match
  position_match <- sapply(seq(ncol(alignment)),
                                  function(x) {
                                    mean(alignment[,x] == reference_segment_string[x])
                                  })
  return(list(percent_match = percent_match, position_match = position_match))
}

for(segment_number in 1:8) {
  segment_name <- paste0("segment", segment_number)
  print(segment_name)
  # compute matches
  assign(segment_name,
         readDNAStringSet(file.path(B_ref_dir, "influenzaB_human_2015-2021",
                                    fasta_files[segment_number])))
  print(paste(segment_name, ":", "alignment to B/Lee/1940"))
  assign(paste0(segment_name, "_lee"),
         calculate_match(get(segment_name), B_Lee[segment_number]))
  print(paste(segment_name, ":", "alignment to consensus"))
  assign(paste0(segment_name, "_consensus"),
         calculate_match(get(segment_name), B_consensus[segment_number]))
  # match across segment
  assign(paste0(segment_name, "_percentMatch_data"),
         data.frame(lee=get(paste0(segment_name, "_lee"))$percent_match,
                    consensus=get(paste0(segment_name, "_consensus"))$percent_match,
                    segment=segment_name))
  # assign(paste0(segment_name, "_percentMatch_plot"),
  #        ggplot(data.frame(lee=get(paste0(segment_name, "_lee"))$percent_match,
  #                          consensus=get(paste0(segment_name, "_consensus"))$percent_match),
  #               aes(x=consensus, y=lee)) +
  #          geom_hex() + geom_abline(slope=1, intercept=0) +
  #          theme_bw() + scale_fill_gradient(low="grey", high="blue") +
  #          xlab("% match to consensus (2015-2021)") +
  #          ylab("% match to B/Lee/1940") +
  #          ggtitle(paste("Influenza B:", segment_name),
  #                  subtitle="% match across sample genome") +
  #          labs(fill="# samples") + xlim(0,1) + ylim(0,1))
  # match across position
  if(width(B_Lee[segment_number]) > width(B_consensus[segment_number])) {
    positions <- pattern(pairwiseAlignment(B_consensus[segment_number],
                                           B_Lee[segment_number]))
    positions <- unlist(strsplit(as.character(positions), split=""))
    positions <- which(positions != "-")
    assign(paste0(segment_name, "_positionMatch_data"),
           data.frame(lee=get(paste0(segment_name, "_lee"))$position_match[positions],
                      consensus=get(paste0(segment_name, "_consensus"))$position_match,
                      segment=segment_name))
    # assign(paste0(segment_name, "_positionMatch_plot"),
    #        ggplot(data.frame(lee=get(paste0(segment_name, "_lee"))$position_match[positions],
    #                          consensus=get(paste0(segment_name, "_consensus"))$position_match),
    #               aes(x=consensus, y=lee)) +
    #          geom_hex() + geom_abline(slope=1, intercept=0) +
    #          theme_bw() + scale_fill_gradient(low="grey", high="blue") +
    #          xlab("% match to consensus (2015-2021)") + ylab("% match to B/Lee/1940") +
    #          ggtitle("", subtitle="% match across position") +
    #          labs(fill="# positions") + xlim(0,1) + ylim(0,1))
  } else {
    positions <- pattern(pairwiseAlignment(B_Lee[segment_number],
                                           B_consensus[segment_number]))
    positions <- unlist(strsplit(as.character(positions), split=""))
    positions <- which(positions != "-")
    assign(paste0(segment_name, "_positionMatch_data"),
           data.frame(lee=get(paste0(segment_name, "_lee"))$position_match,
                      consensus=get(paste0(segment_name, "_consensus"))$position_match[positions],
                      segment=segment_name))
    # assign(paste0(segment_name, "_positionMatch_plot"),
    #        ggplot(data.frame(lee=get(paste0(segment_name, "_lee"))$position_match,
    #                          consensus=get(paste0(segment_name, "_consensus"))$position_match[positions]),
    #               aes(x=consensus, y=lee)) +
    #          geom_hex() + geom_abline(slope=1, intercept=0) +
    #          theme_bw() + scale_fill_gradient(low="grey", high="blue") +
    #          xlab("% match to consensus (2015-2021)") + ylab("% match to B/Lee/1940") +
    #          ggtitle("", subtitle="% match across position") +
    #          labs(fill="# positions") + xlim(0,1) + ylim(0,1))
  }
  # assign(paste0(segment_name, "_plot"),
  #        get(paste0(segment_name, "_percentMatch_plot")) +
  #          get(paste0(segment_name, "_positionMatch_plot")))
}

all_position_match <- lapply(grep("positionMatch", ls(), value=T), get)
all_position_match <- do.call(rbind, all_position_match)

all_percent_match <- lapply(grep("percentMatch", ls(), value=T), get)
all_percent_match <- do.call(rbind, all_percent_match)

summary_plot <- (ggplot(all_percent_match, aes(x=consensus, y=lee)) +
                   geom_hex() + geom_abline(slope=1, intercept=0) + facet_grid(~segment) +
                   theme_bw() + scale_fill_gradient(low="lightgrey", high="blue") +
                   xlab("% match to consensus (2015-2021)") + ylab("% match to B/Lee/1940") +
                   ggtitle("", subtitle="% match across sample genome")) /
  (ggplot(all_position_match, aes(x=consensus, y=lee)) +
     geom_hex() + geom_abline(slope=1, intercept=0) + facet_grid(~segment) +
     theme_bw() + scale_fill_gradient(low="lightgrey", high="blue") +
     xlab("% match to consensus (2015-2021)") + ylab("% match to B/Lee/1940") +
     ggtitle("", subtitle="% match across positions"))
ggsave(file.path(here(), "outputs", "influenzaB_human_2015-2020_cas13a_neg_20nt",
                 "plot_consensus_v_lee1940.pdf"), summary_plot,
       device="pdf", width=10, height=5, units="in")
