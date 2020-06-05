##################################################
### GUI to visualize guide screening

library(shiny)
library(ggplot2)
library(here)
library(patchwork)

# import data -------------------------------------------------------------

plusStrand <- read.table(file.path(here(), "outputs/cas13a_20nt/cas13a_results_summary.txt"))
minusStrand <- read.table(file.path(here(), "outputs/cas13a_minusStrand_20nt/cas13a_minusStrand_results_summary.txt"))
dat <- rbind(plusStrand, minusStrand)

# filter candidate guides
filtered_dat <- subset(dat, has_crRNA_hairpin) # n = 33008 (+: 17070; -: 115938)
filtered_dat <- subset(filtered_dat, crRNA_spacer_basepairs <= 4) # n = 17229 (+: 9159; -: 8070)
filtered_dat <- subset(filtered_dat, specificity==1) # n = 11721 (+: 6219, -: 5502)
filtered_dat <- subset(filtered_dat, match_against_hg38==0) # n = 7104 (+: 3729; -: 3375)
filtered_dat <- subset(filtered_dat, antitag != "GUUU") # n = 7047 (+: 3702; -: 3345)
filtered_dat <- subset(filtered_dat, sensitivity_01 >= 0.95) # n = 6969 (+: 3678; -: 3291)

# load viral abundance coverage data
coverage_id <- "PRJNA616446"
abundance_binSize <- 300
id_cov <- read.table(file.path(here(), "ref_data/RNA_expression", paste0(coverage_id, "_mapped.cov")),
                     header=F, stringsAsFactors=F, col.names=c("genome", "pos", "coverage"))
id_cov$bin <- cut(id_cov$pos, breaks=ceiling(max(id_cov$pos)/abundance_binSize), labels=F)
id_cov <- aggregate(coverage~bin, data=id_cov, FUN=median)
id_cov$bin <- (id_cov$bin-1)*abundance_binSize + abundance_binSize/2

# filter 
filters <- c("All guides" = "dat",
             "Filtered guides" = "filtered_dat",
             "Ordered guides" = "ordered_dat")
filter_label <- names(filters)
names(filter_label) <- filters

# features
features <- c("Sensitivity (no mismatch)" = "sensitivity_0",
              "Sensitivity (1 mismatch)" = "sensitivity_01",
              "Specificity" = "specificity",
              "GC content" = "GC_content",
              "Spacer base-pairing" = "crRNA_spacer_basepairs",
              "crRNA hairpin intact" = "has_crRNA_hairpin",
              "Alignments to human transcriptome" = "match_against_hg38")
feature_label <- names(features)
names(feature_label) <- features
feature_label["match_against_hg38"] <- "# alignments to \nhuman transcriptome"
feature_desc <- c("% SARS-CoV-2 samples detected (no mismatches)",
                  "% SARS-CoV-2 samples detected (1 mismatch)",
                  "% human coronaviruses missed (2 mismatches)",
                  "Spacer GC content",
                  "# basepairs in spacer",
                  "crRNA hairpin intact",
                  "# alignments to human transcriptome")
names(feature_desc) <- features

# columns
display_columns <- c("start", "strand", "target", "spacer", "antitag", "sensitivity_01", "specificity")

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
    return(datatable(plot_data[, display_columns], rownames=F, selection="multiple"))
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
      scale_color_gradient2(high="red", low="blue", mid="grey", midpoint=mean(plot_data[,input$var3]),
                            name=feature_label[input$var3]) +
      ggtitle("Guide design scores", subtitle=paste("n =", nrow(plot_data)))
    # univariate plots
    var1_plot <- ggplot(plot_data, aes_string(input$var1)) + geom_density(kernel="gaussian", fill=2, col=2) +
      theme_bw() + ggtitle(feature_label[input$var1]) + xlab(feature_desc[input$var1])
    var2_plot <- ggplot(plot_data, aes_string(input$var2)) + geom_density(kernel="gaussian", fill=3, col=3) +
      theme_bw() + ggtitle(feature_label[input$var2]) + xlab(feature_desc[input$var2])
    var3_plot <- ggplot(plot_data, aes_string(input$var3)) + geom_density(kernel="gaussian", fill=4, col=4) +
      theme_bw() + ggtitle(feature_label[input$var3]) + xlab(feature_desc[input$var3])
    return(scatter_plot + (var1_plot / var2_plot / var3_plot) + plot_layout(widths=c(2,1)))
  })
  
}


# run app -----------------------------------------------------------------

shinyApp(ui=ui, server=server)
