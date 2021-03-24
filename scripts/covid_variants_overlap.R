### check for overlap between NCR guides and COVID variants
# variants: B.1.1.7, B.1.351, P.1

library(here)

variants <- c(3268, 5389, 6955, 11289, 21766, 21992, 23064, 23272, 23605, 24507, 
              24915, 27973, 28049, 28112, 28281, 28978, 23403, 21801, 22206, 
              23013, 23064, 23664, 21615, 22299, 22812, 22287, 241, 100, 28882, 
              28883, 28884, 29149, 28976, 28629, 29754, 3038, 14410, 12054, 
              11825, 10668, 12965, 26621, 28254, 23404, 25089, 23013, 22029,
              22779, 22920, 22878, 22974)

plate1 <- read.table(file.path(here(), "outputs", "ordered_guides.bed"), 
                     col.names=c("chr", "start", "end", "id"))
plate2 <- read.csv(file.path(here(), "outputs", "control_guides_20210106.csv"))
plate2$id <- paste0("NCR_", 1305:1400)

overlap <- cbind(variant_pos = variants,
                 plate1 = sapply(variants, 
                                 function(x) {
                                   paste(subset(plate1, 
                                                (start <= x) & ((start + 19) >= x))$id, 
                                         collapse=", ")
                                 }),
                 plaste2 = sapply(variants, 
                                  function(x) {
                                    paste(subset(plate2, 
                                                 (start <= x) & ((start + 19) >= x))$id, 
                                          collapse=", ")
                                  }))
