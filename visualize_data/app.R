##################################################
### GUI to visualize guide screening

library(shiny)
library(ggplot2)
library(patchwork)
# library(ggbio)
# library(GenomicRanges)

# import data -------------------------------------------------------------

plusStrand <- read.table("data/cas13a_results_summary.txt")
minusStrand <- read.table("data/cas13a_minusStrand_results_summary.txt")
dat <- rbind(plusStrand, minusStrand)

offtarget_candidates <- read.table("data/candidate_spacers.txt", stringsAsFactors=F, header=T)
offtarget_vir <- read.csv("data/20200615_1056_guides.csv", stringsAsFactors=F)
offtarget_bac <- read.table("data/mat_refbac_20200609.tsv", header=T, stringsAsFactors=F)
offtarget_bac$hits_01 <- offtarget_bac$X0 + offtarget_bac$X1

amplicons <- readLines("data/winner_amplicons_mapped.sam")
amplicons <- amplicons[!grepl("^@", amplicons)]
amplicons <- data.frame(matrix(unlist(strsplit(amplicons, split="\t")), nrow=length(amplicons), byrow=T), stringsAsFactors=F)
amplicons <- amplicons[,c(1, 4, 10)]
colnames(amplicons) <- c("amplicon", "pos", "seq")
amplicons$pos <- as.numeric(as.character(amplicons$pos))
amplicons$end <- amplicons$pos + nchar(amplicons$seq)

# add offtarget data
dat$offtarget_vir <- NA
dat$offtarget_vir[prodlim::row.match(offtarget_candidates[, c("start", "strand")], dat[, c("start", "strand")])] <- T
dat$offtarget_vir[prodlim::row.match(offtarget_vir[, c("start", "strand")], dat[, c("start", "strand")])] <- F
dat$offtarget_bac <- NA
dat$offtarget_bac[prodlim::row.match(offtarget_bac[, c("start_position", "strand")], 
                                     dat[, c("start", "strand")])] <- (offtarget_bac$hits_01>0)
dat$offtarget_any <- !(!dat$offtarget_bac & !dat$offtarget_vir)
offtarget_dat <- subset(dat, !offtarget_any)

# # convert logicals to numerics
# dat$has_crRNA_hairpin <- as.numeric(dat$has_crRNA_hairpin)
# dat$offtarget_any <- as.numeric(dat$offtarget_any)
# dat$offtarget_bac <- as.numeric(dat$offtarget_bac)
# dat$offtarget_vir <- as.numeric(dat$offtarget_vir)

# filter candidate guides: QC
qc_dat <- subset(dat, has_crRNA_hairpin==1) # n = 33008 (+: 17070; -: 115938)
qc_dat <- subset(qc_dat, crRNA_spacer_basepairs <= 4) # n = 17229 (+: 9159; -: 8070)
qc_dat <- subset(qc_dat, specificity==1) # n = 11721 (+: 6219, -: 5502)
qc_dat <- subset(qc_dat, match_against_hg38==0) # n = 7104 (+: 3729; -: 3375)
qc_dat <- subset(qc_dat, antitag != "GUUU") # n = 7047 (+: 3702; -: 3345)
qc_dat <- subset(qc_dat, sensitivity_01 >= 0.95) # n = 6969 (+: 3678; -: 3291)

# # load viral abundance coverage data
# coverage_id <- "PRJNA616446"
# abundance_binSize <- 300
# id_cov <- read.table(file.path(here(), "ref_data/RNA_expression", paste0(coverage_id, "_mapped.cov")),
#                      header=F, stringsAsFactors=F, col.names=c("genome", "pos", "coverage"))
# id_cov$bin <- cut(id_cov$pos, breaks=ceiling(max(id_cov$pos)/abundance_binSize), labels=F)
# id_cov <- aggregate(coverage~bin, data=id_cov, FUN=median)
# id_cov$bin <- (id_cov$bin-1)*abundance_binSize + abundance_binSize/2

# filter 
filters <- c("All guides" = "dat",
             "QC-filtered guides (n=6969)" = "qc_dat",
             "Offtarget-filtered guides (n=258)" = "offtarget_dat")
filter_label <- names(filters)
names(filter_label) <- filters

# features
features <- c("Sensitivity (no mismatch)" = "sensitivity_0",
              "Sensitivity (1 mismatch)" = "sensitivity_01",
              "Specificity" = "specificity",
              "GC content" = "GC_content",
              "Spacer base-pairing" = "crRNA_spacer_basepairs",
              "crRNA hairpin intact" = "has_crRNA_hairpin",
              "Alignments to human transcriptome" = "match_against_hg38",
              "Offtarget (all)" = "offtarget_any",
              "Offtarget (viral)" = "offtarget_vir",
              "Offtarget (bacterial)" = "offtarget_bac")
feature_label <- names(features)
names(feature_label) <- features
feature_label["match_against_hg38"] <- "# alignments to \nhuman transcriptome"
feature_desc <- c("% SARS-CoV-2 samples detected (no mismatches)",
                  "% SARS-CoV-2 samples detected (1 mismatch)",
                  "% human coronaviruses missed (2 mismatches)",
                  "Spacer GC content",
                  "# basepairs in spacer",
                  "crRNA hairpin intact",
                  "# alignments to human transcriptome",
                  "Guide hits any offtarget (1 mismatch)",
                  "Guide hits viral offtarget (1 mismatch)",
                  "Guide hits bacterial offtarget (1 mismatch)")
names(feature_desc) <- features
logical_features <- colnames(dat)[sapply(colnames(dat), function(x) class(dat[,x])=="logical")]

# columns
display_columns <- c("start", "strand", "target", "spacer", "antitag", "sensitivity_01", "offtarget_any")

# define UI ---------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("COVID-19 Guide Design & Evaluation"),
  
  sidebarLayout(
    sidebarPanel(h2("Select data"),
                 selectInput("data_filter", "Select data:", filters),
                 radioButtons("strand_filter", "Strand Filter", c("both strands" = "both",
                                                                  "+ strand" = "plus",
                                                                  "- strand" = "minus")),
                 
                 h2("Select features"),
                 selectInput("var1", "Feature 1:", features),
                 selectInput("var2", "Feature 2:", features),
                 selectInput("var3", "Feature 3:", features),
                 
                 h2("Select plot region"),
                 sliderInput("genome_range", "Genomic range",
                             min=1, max=nrow(plusStrand), value=c(1, nrow(plusStrand)))
    ),
    
    mainPanel(
      # plotOutput("genome_plot"),
      plotOutput("scatter_plot"),
      p(class="text-center", downloadButton('x3', "Download Data"),
        DT::dataTableOutput("print_table"))
    )
  )
)

# define server logic -----------------------------------------------------

server <- function(input, output) {
  
  output$print_table <- DT::renderDataTable({
    plot_data <- get(input$data_filter)
    plot_data <- subset(plot_data, start >= input$genome_range[1])
    plot_data <- subset(plot_data, start <= input$genome_range[2]-20)
    plot_data <- switch(input$strand_filter,
                        "both" = plot_data,
                        "plus" = subset(plot_data, strand=="+"),
                        "minus" = subset(plot_data), strand=="-")
    return(DT::datatable(plot_data[, display_columns], rownames=F, selection="multiple"))
  })
  
  output$scatter_plot <- renderPlot({
    # subset data
    plot_data <- get(input$data_filter)
    plot_data <- subset(plot_data, start >= input$genome_range[1])
    plot_data <- subset(plot_data, start <= input$genome_range[2]-20)
    plot_data <- switch(input$strand_filter,
                        "both" = plot_data,
                        "plus" = subset(plot_data, strand=="+"),
                        "minus" = subset(plot_data), strand=="-")
    # scatter plot
    scatter_plot <- ggplot(plot_data, aes_string(input$var1, input$var2, col=input$var3)) +
      geom_jitter(width=(max(plot_data[input$var1])-min(plot_data[input$var1]))/20,
                  height=(max(plot_data[input$var1])-min(plot_data[input$var1]))/20,
                  alpha=0.5) +
      theme_bw() + xlab(feature_label[input$var1]) + ylab(feature_label[input$var2]) +
      ggtitle("Guide design scores", subtitle=paste("n =", nrow(plot_data)))
    if(!(input$var3 %in% logical_features)) {
      scatter_plot <- scatter_plot + scale_color_gradient2(high="red", low="blue", mid="grey", 
                                                           midpoint=mean(plot_data[,input$var3]),
                                                           name=feature_label[input$var3])
    }
    # univariate plots
    var1_plot <- ggplot(plot_data, aes_string(input$var1)) + theme_bw() + 
      ggtitle(feature_label[input$var1]) + xlab(feature_desc[input$var1])
    if(input$var1 %in% logical_features) { var1_plot <- var1_plot + geom_bar() } 
    else { var1_plot <- var1_plot + geom_density(kernel="gaussian", fill=2, col=2) }
    var2_plot <- ggplot(plot_data, aes_string(input$var2)) + theme_bw() + 
      ggtitle(feature_label[input$var2]) + xlab(feature_desc[input$var2])
    if(input$var2 %in% logical_features) { var2_plot <- var2_plot + geom_bar() } 
    else { var2_plot <- var2_plot + geom_density(kernel="gaussian", fill=3, col=3) }
    var3_plot <- ggplot(plot_data, aes_string(input$var3)) + theme_bw() + 
      ggtitle(feature_label[input$var3]) + xlab(feature_desc[input$var3])
    if(input$var3 %in% logical_features) { var3_plot <- var3_plot + geom_bar() } 
    else { var3_plot <- var3_plot + geom_density(kernel="gaussian", fill=4, col=4) }
    return(scatter_plot + (var1_plot / var2_plot / var3_plot) + plot_layout(widths=c(2,1)))
  })
  
  # output$genome_plot <- renderPlot({
  #   # subset data
  #   plot_data <- get(input$data_filter)
  #   plot_data <- subset(plot_data, start >= input$genome_range[1])
  #   plot_data <- subset(plot_data, start <= input$genome_range[2]-20)
  #   plot_data <- switch(input$strand_filter,
  #                       "both" = plot_data,
  #                       "plus" = subset(plot_data, strand=="+"),
  #                       "minus" = subset(plot_data), strand=="-")
  #   plus_strand <- subset(plot_data, strand=="+")
  #   minus_strand <- subset(plot_data, strand=="-")
  #   amplicons_subset <- subset(amplicons, pos >= input$genome_range[1] & end <= input$genome_range[2])
  #   # make GRanges objects
  #   if(nrow(plus_strand) > 0) {
  #     plus_strand_obj <- GRanges(seqnames="SARS-CoV-2",
  #                                ranges=IRanges(start=plus_strand$start, width=20, 
  #                                               names=paste("+", plus_strand$start, sep="_")),
  #                                type="target", strand="+")
  #   }
  #   if(nrow(minus_strand) > 0) {
  #     minus_strand_obj <- GRanges(seqnames="SARS-CoV-2",
  #                                 ranges=IRanges(start=minus_strand$start, width=20, 
  #                                                names=paste("-", minus_strand$start, sep="_")),
  #                                 type="target", strand="-")
  #   }
  #   if(nrow(amplicons_subset) > 0) {
  #     amplicons_obj <- GRanges(seqnames="SARS-CoV-2", 
  #                              ranges=IRanges(start=amplicons_subset$pos, end=amplicons_subset$end,
  #                                             name=amplicons_subset$amplicon),
  #                              type="amplicon")
  #   }
  #   # generate plot
  #   if(nrow(amplicons_subset) > 0) {
  #     if(nrow(plus_strand) > 0 & nrow(minus_strand) > 0) {
  #       genome_plot <- tracks("+ strand" = autoplot(plus_strand_obj, fill=2),
  #                             "- strand" = autoplot(minus_strand_obj, fill=4),
  #                             amplicon = autoplot(amplicons_obj),
  #                             heights=c(1,1,0.2))
  #     } else if(nrow(plus_strand) > 0 & nrow(minus_strand) == 0) {
  #       genome_plot <- tracks("+ strand" = autoplot(plus_strand_obj, fill=2),
  #                             amplicon = autoplot(amplicons_obj),
  #                             heights=c(1,0.2))
  #     } else if(nrow(plus_strand) == 0 & nrow(minus_strand) > 0) {
  #       genome_plot <- tracks("- strand" = autoplot(minus_strand_obj, fill=4),
  #                             amplicon = autoplot(amplicons_obj),
  #                             heights=c(1,0.2))
  #     } else {
  #       genome_plot <- tracks(amplicon = autoplot(amplicons_obj))
  #     }
  #   } else {
  #     if(nrow(plus_strand) > 0 & nrow(minus_strand) > 0) {
  #       genome_plot <- tracks("+ strand" = autoplot(plus_strand_obj, fill=2),
  #                             "- strand" = autoplot(minus_strand_obj, fill=4))
  #     } else if(nrow(plus_strand) > 0 & nrow(minus_strand) == 0) {
  #       genome_plot <- tracks("+ strand" = autoplot(plus_strand_obj, fill=2))
  #     } else if(nrow(plus_strand) == 0 & nrow(minus_strand) > 0) {
  #       genome_plot <- tracks("- strand" = autoplot(minus_strand_obj, fill=4))
  #     } else {
  #       return(Ideogram(genome="wuhCor1"))
  #     }
  #   }
  #   return(genome_plot + theme_bw() + xlim(input$genome_range[1], input$genome_range[2]))
  # })
  
}


# run app -----------------------------------------------------------------

shinyApp(ui=ui, server=server)
