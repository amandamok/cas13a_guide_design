##################################################
### generate reference influenza genome segments

library(here)
library(Biostrings)
library(msa)

ref_dir <- file.path(here(), "ref_data")
window_size <- 20

output_dir <- file.path(here(), "outputs", "flu_cas13a_neg_20nt")
if(!dir.exists(output_dir)) { dir.create(output_dir) }

a_california_strain <- "A/California/04/2009"
b_brisbane_strain <- "B/Birsbane/60/2008"
b_victoria_strain <- "B/Victoria/504/2000"

a_california_file <- gsub("\\/", "_", a_california_strain)
b_brisbane_file <- gsub("\\/", "_", b_brisbane_strain)
b_victoria_file <- gsub("\\/", "_", b_victoria_strain)

# define functions --------------------------------------------------------

get_annotation <- function(fasta_headers) {
  # fasta_headers: character vector from .fasta file
  annotations <- sapply(fasta_headers,
                        function(x) {
                          x <- strsplit(x, split="\\|")[[1]]
                          sapply(c("gb", "Strain Name", "Segment"),
                                 function(key) {
                                   sub(paste0(key, ":"), "",
                                       grep(key, x, value=T))
                                 })
                        })
  annotations <- data.frame(t(annotations),
                            row.names=NULL, stringsAsFactors=F)
  colnames(annotations) <- c("name", "strain", "segment")
  annotations$segment <- as.numeric(annotations$segment)
  return(annotations)
}

trim_nt <- function(consensus_string) {
  # consensus_string: character
  consensus_string <- sub("^*(\\?|-)+", "", consensus_string)
  consensus_string <- sub("(\\?|-)+$", "", consensus_string)
  return(consensus_string)
}

get_consensus <- function(segments, annotations) {
  # segments: DNAStringSet; genome segment sequences
  # annotations: data.frame; output from get_annotation()
  sapply(1:8,
         function(x) {
           segment <- segments[annotations$segment==x]
           gsub("\\?", "N", trim_nt(msaConsensusSequence(msa(segment))))
         })
}

get_targets <- function(segments, strain_name, window_width=20) {
  # segments: character vector; consensus sequences from get_consensus()
  # window_width: integer; desired length of Cas13a spacer
  # strain_name: character; name of strain
  windows <- lapply(1:8,
                    function(x) {
                      segment <- segments[x]
                      num_windows <- nchar(segment) - window_width - 4 + 1 # 4nt for antitag
                      targets <- sapply(seq(num_windows),
                                        function(y) {
                                          substr(segment,
                                                 y, y + window_width - 1)
                                        })
                      spacers <- reverse_complement(targets)
                      antitags <-sapply(seq(num_windows),
                                        function(y) {
                                          substr(segment,
                                                 y + window_width,
                                                 y + window_width + 4 - 1)
                                        })
                      return(data.frame(segment = x,
                                        position = 1:num_windows,
                                        target = targets,
                                        spacer = spacers,
                                        antitag = antitags))
                    })
  windows <- data.frame(strain = strain_name, do.call(rbind, windows))
  return(windows)
}

reverse_complement <- function(sequences) {
  # sequences: character vector; DNA strings to be converted to RNA strings
  sequences <- Biostrings::DNAStringSet(sequences)
  sequences <- Biostrings::reverseComplement(sequences)
  sequences <- as.character(sequences)
  sequences <- gsub("T", "U", sequences)
  return(sequences)
}

# load genome segments ----------------------------------------------------

a_california_fname <- file.path(ref_dir,
                                paste0("flu_", a_california_file, ".fasta"))
a_california_seq <- readDNAStringSet(a_california_fname)
a_california_anno <- get_annotation(names(a_california_seq))

b_brisbane_fname <- file.path(ref_dir,
                              paste0("flu_", b_brisbane_file, ".fasta"))
b_brisbane_seq <- readDNAStringSet(b_brisbane_fname)
b_brisbane_anno <- get_annotation(names(b_brisbane_seq))

b_victoria_fname <- file.path(ref_dir,
                              paste0("flu_", b_victoria_file, ".fasta"))
b_victoria_seq <- readDNAStringSet(b_victoria_fname)
b_victoria_anno <- get_annotation(names(b_victoria_seq))

# align segments ----------------------------------------------------------

a_california_segments <- get_consensus(a_california_seq, a_california_anno)
b_brisbane_segments <- get_consensus(b_brisbane_seq, b_brisbane_anno)
b_victoria_segments <- get_consensus(b_victoria_seq, b_victoria_anno)

# break segments into 20nt windows ----------------------------------------

a_california <- get_targets(a_california_segments, a_california_strain)
b_brisbane <- get_targets(b_brisbane_segments, b_brisbane_strain)
b_victoria <- get_targets(b_victoria_segments, b_victoria_strain)

# write outputs -----------------------------------------------------------

windows <- rbind(a_california, b_brisbane, b_victoria)
windows$GC_content <- Biostrings::letterFrequency(Biostrings::RNAStringSet(windows$spacer),
                                                  letters="GC") / window_size
windows$U_content <- Biostrings::letterFrequency(Biostrings::RNAStringSet(windows$spacer),
                                                 letters="U") / window_size

# write output table
write.table(windows, file=file.path(output_dir, "windows.txt"),
            quote=F, sep="\t", row.names=F)

# write target sequences
writeLines(as.character(windows$target), con=file.path(output_dir, "targets.txt"))

# write spacer sequences
writeLines(as.character(windows$spacer), con=file.path(output_dir, "spacers.txt"))

# write target sequences as fasta file
windows_fasta <- DNAStringSet(windows$target)
names(windows_fasta) <- with(windows, paste(strain, segment, position, sep="_"))
writeXStringSet(windows_fasta, filepath=file.path(output_dir, "targets.fa"))

# write strain consensus sequences
a_california_segments <- DNAStringSet(a_california_segments)
names(a_california_segments) <- paste(a_california_strain, seq_along(a_california_segments), sep="_")
writeXStringSet(a_california_segments, filepath=file.path(output_dir,
                                                          paste0(a_california_file, ".fasta")))
b_brisbane_segments <- DNAStringSet(b_brisbane_segments)
names(b_brisbane_segments) <- paste(b_brisbane_strain, seq_along(b_brisbane_segments), sep="_")
writeXStringSet(b_brisbane_segments, filepath=file.path(output_dir,
                                                          paste0(b_brisbane_file, ".fasta")))
b_victoria_segments <- DNAStringSet(b_victoria_segments)
names(b_victoria_segments) <- paste(b_victoria_strain, seq_along(b_victoria_segments), sep="_")
writeXStringSet(b_victoria_segments, filepath=file.path(output_dir,
                                                          paste0(b_victoria_file, ".fasta")))

q(save="no")