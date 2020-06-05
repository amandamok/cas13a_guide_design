##################################################
### GUI to visualize guide screening

library(shiny)
library(ggplot2)
library(here)

source(file.path(here(), "scripts/helper.R"))

coverage_id <- "PRJNA616446"
abundance_binSize <- 300
abundance <- load_coverage(coverage_id, abundance_binSize)

# import data -------------------------------------------------------------

plusStrand <- read.table(file.path(here(), "outputs/cas13a_20nt/cas13a_results_summary.txt"), stringsAsFactors=F)
minusStrand <- read.table(file.path(here(), "outputs/cas13a_minusStrand_20nt/cas13a_minusStrand_results_summary.txt"), stringsAsFactors=F)
dat <- rbind(plusStrand, minusStrand)

# filter candidate guides
filtered_dat <- subset(dat, has_crRNA_hairpin) # n = 33008 (+: 17070; -: 115938)
filtered_dat <- subset(filtered_dat, crRNA_spacer_basepairs <= 4) # n = 17229 (+: 9159; -: 8070)
filtered_dat <- subset(filtered_dat, specificity==1) # n = 11721 (+: 6219, -: 5502)
filtered_dat <- subset(filtered_dat, match_against_hg38==0) # n = 7104 (+: 3729; -: 3375)
filtered_dat <- subset(filtered_dat, antitag != "GUUU") # n = 7047 (+: 3702; -: 3345)
filtered_dat <- subset(filtered_dat, sensitivity_01 >= 0.95) # n = 6969 (+: 3678; -: 3291)

# define UI ---------------------------------------------------------------

filters <- c("All guides" = "dat",
             "Filtered guides" = "filtered_dat",
             "Ordered guides" = "ordered_dat")
features <- c("Sensitivity (no mismatch)" = "sensitivity_0",
              "Sensitivity (1 mismatch)" = "sensitivity_01",
              "Specificity" = "specificity",
              "GC content" = "GC_content",
              "Spacer base-pairing" = "crRNA_spacer_basepairs",
              "crRNA hairpin intact" = "has_crRNA_hairpin",
              "Alignments to human transcriptome" = "match_against_hg38")

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
      # img(src="rstudio.png", height=140, width=400),
      # plotOutput("diagnostic_plot"),
      textOutput("selected_var")
    )
  )
)

# define server logic -----------------------------------------------------

filter_label <- names(filters)
names(filter_label) <- filters
feature_label <- names(features)
names(feature_label) <- features
feature_desc <- c("% SARS-CoV-2 samples detected (no mismatches)",
                  "% SARS-CoV-2 samples detected (1 mismatch)",
                  "% human coronaviruses missed (2 mismatches)",
                  "Spacer GC content",
                  "# basepairs in spacer",
                  "crRNA hairpin intact",
                  "# alignments to human transcriptome")
names(feature_desc) <- features

server <- function(input, output) {
  
  # output$diagnostic_plot <- renderPlot(
  #   plot_diagnostic(get(input$data_filter), abundance, filter=filter_label[input$data_filter], 
  #                   alpha=0.5, jitter_x=0.001, jitter_y=0.05, 
  #                   var1=input$var1, var1_name=feature_label[input$var1], var1_desc=feature_desc[input$var1],
  #                   var2=input$var2, var2_name=feature_label[input$var2], var2_desc=feature_desc[input$var2],
  #                   var3=input$var3, var3_name=feature_label[input$var3], var3_desc=feature_desc[input$var3]))

  output$selected_var <- renderText({"You have selected this"})
}


# run app -----------------------------------------------------------------

shinyApp(ui=ui, server=server)
