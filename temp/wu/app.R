# Reference: https://www.nextflow.io/docs/latest/getstarted.html
## Shinyapps.io detects and installs packages for you automatically when you call deployApp(). 
## Do not need, nor should have any calls to install.packages() as below anywhere in your source code.
## Below installation check is only for local installation
# # 1. CRAN packages----
# packages <- c("shinydashboard",
#               "data.table",
#               "dplyr",
#               "DT",
#               "farver",
#               "fgsea",
#               "ggdendro",
#               "ggplot2",
#               "gridExtra",
#               "knitr",
#               "MASS",
#               "packrat",
#               "pheatmap",
#               "plotly",
#               "RColorBrewer",
#               "readxl",
#               "shiny",
#               "shinyFiles",
#               "shinythemes",
#               "shinyWidgets",
#               "stringr",
#               "tibble",
#               "units",
#               "VennDiagram",
#               "zip",
#               "tidyverse")
# if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(packages, rownames(installed.packages())))
# }
#
# # 2. Bioconductor packages----
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocInstaller", version = "3.8")
# BiocManager::install("DESeq2", version = "3.8")
# BiocManager::install("DEGseq", version = "3.8")
# BiocManager::install("GOSemSim", version = "3.8")
# BiocManager::install("ChIPseeker", version = "3.8")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", version = "3.8")
# BiocManager::install("DSS", version = "3.8")
# BiocManager::install("farver", version = "3.8")
# BiocManager::install("units", version = "3.8")
# BiocManager::install("fgsea", version = "3.8")
# BiocManager::install("org.Mm.eg.db", version = "3.8")

options(stringsAsFactors = FALSE)
options(shiny.maxRequestSize = 30*1024^2) 
options(repos = BiocInstaller::biocinstallRepos())

library(shinydashboard)
library(DT)
library(shiny)
library(shinythemes)
library(shinyFiles)
library(shinyWidgets)



library(data.table)
library(ggdendro)
library(ggplot2)
library(gridExtra)
library(knitr)
library(MASS)
library(pheatmap)
library(plotly)
library(RColorBrewer)
library(readxl)
library(stringr)
library(tibble)
library(VennDiagram)
library(zip)

library(BiocInstaller)
library(DESeq2)
library(DEGseq)
library(BiocParallel)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DSS)

# load dplyr last to avoid "select" being masked from biomaRt or plotly
library(tidyverse)

# Determine the OS
sysOS <- Sys.info()[['sysname']]

# Specify root folder----
if (sysOS == "Linux") {
  volumes <- c("Home" = "/home/administrator/",
               "Data Storage" = "/datastorage/")
} else if (sysOS == "Windows") {
  volumes <- c("Local Disk (C:)" = "C:/")
}

# source functions
source("source/plots.R")
source("source/data_prep.R")

ui <- dashboardPage(dashboardHeader(title = "NGS Pipeline"),
                    dashboardSidebar(sidebarMenu(menuItem(text = "Introduction", 
                                                          tabName = "introduction", 
                                                          icon = icon("dashboard")),
                                                 menuItem(text = "FASTQ Processing", 
                                                          tabName = "count", 
                                                          icon = icon("th")),
                                                 menuItem(text = "RNA-seq Analysis", 
                                                          tabName = "rna-seq_analysis", 
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA MethylSeq Analysis",
                                                          tabName = "dna-seq_analysis",
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA vs RNA",
                                                          tabName = "dna_vs_rna",
                                                          icon = icon("th")))),
                    dashboardBody(tabItems(
                      
                      ## ---------- Introduction ------------
                      tabItem(tabName = "introduction",
                              h2("Instructions"),
                              p("This interactive web application (NGS pipeline) is developed in R with Shiny to 
                                1) Align RNA and DNA-seq
                                2) conduct differential expression (DE) analysis with ",
                                a("DEGseq, ",
                                  href = "https://bioconductor.org/packages/release/bioc/html/DEGseq.html"),
                                a("DESeq2 ",
                                  href = "http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
                                " based on the provided input data (raw count table and experimental design table).",
                                style = "padding-left: 0em"),
                              h3("1. RNA Analysis Data input", 
                                 style = "padding-left: 1em"),
                              h4("1.1 Experimental Design Table",
                                 style = "padding-left: 3em"), 
                              p("This input data contains summarized experimental design information for each sample. 
                                The first column, Sample Name, will be the same as in Count Table. The Second column, 
                                Sample Label, will be shown in the result plots and tables. 
                                A template of the design table can be downloaded ",
                                a("here.", 
                                  href = "https://drive.google.com/uc?export=download&id=1iWeDU8JK5mZxpp6fd9ipf9JXuKKztv7Y")
                                , style = "padding-left: 5em")),
                    
                    ## -------- FASTQ Processing -------------------
                    tabItem(tabName = "count",
                            sidebarPanel(
                              wellPanel(shinyDirButton(id = "directory", 
                                                       label = "Select FASTQ Folder", 
                                                       title = "Please select a folder"),
                                        radioButtons(inputId = "is_rna",
                                                     label = "",
                                                     inline = TRUE,
                                                     choices = list("RNA" = TRUE, 
                                                                    "DNA" = FALSE),
                                                     selected = TRUE),
                                        radioButtons(inputId = "isSingleEnd",
                                                     label = "Sequencing:",
                                                     inline = TRUE,
                                                     choices = list("Single Ended" = "true", 
                                                                    "Pair Ended" = "false"),
                                                     selected = "true"),
                                        
                                        fluidRow(
                                          column(4,
                                                 actionButton(inputId = "align",
                                                              label = "Align"))),
                                        
                                        fluidRow(
                                          checkboxInput(inputId = "show_advanced_para",
                                                        label = "Show Advanced Setting",
                                                        value = FALSE)),
                                        conditionalPanel(condition = "input.show_advanced_para == true",
                                                         textInput(inputId = "job_name",
                                                                   label = "job name (mmddy + x job)",
                                                                   value = paste0(substr(format(Sys.Date(),"%m%d%y"),1,4),
                                                                                  substr(format(Sys.Date(),"%m%d%y"),6,6),
                                                                                  "1")),
                                                         selectInput(inputId = "ntasks",
                                                                     label = "ntasks",
                                                                     choices = c(1),
                                                                     multiple = FALSE),
                                                         selectInput(inputId = "cpus_per_task",
                                                                     label = "cpus per task",
                                                                     choices = c(1:32),
                                                                     multiple = FALSE,
                                                                     selected = 2),
                                                         selectInput(inputId = "mem",
                                                                     label = "RAM",
                                                                     choices = c(paste0(c(16, 32, 64, 120, 192), "GB")),
                                                                     multiple = FALSE,
                                                                     selected = "64GB"),
                                                         selectInput(inputId = "time",
                                                                     label = "total run time limit (HH:MM:SS)",
                                                                     choices = c(paste0(c(2, 4, 6, 8, 12, 24),":00:00")),
                                                                     multiple = FALSE,
                                                                     selected = "2:00:00"),
                                                         textInput(inputId = "output",
                                                                   label = "STDOUT output file",
                                                                   value = "slurm.%N.%j.out"),
                                                         textInput(inputId = "error",
                                                                   label = "STDERR output file (optional)",
                                                                   value = "slurm.%N.%j.err"),
                                                         selectInput(inputId = "partition",
                                                                     label = "partition",
                                                                     choices = c("p_kongt_1", "main"),
                                                                     multiple = FALSE,
                                                                     selected = "p_kongt_1")))),
                            mainPanel(tags$h3("File Directory"),
                                      verbatimTextOutput(outputId = "fileDir"),
                                      tags$h3("Program"),
                                      verbatimTextOutput("program"))
                            # tabItem "count" ends
                            ),
                    
                    ## ---------------- RNA seq analysis -----------------
                    tabItem(tabName = "rna-seq_analysis",
                            tabsetPanel(
                              ## ------- read in data -------------
                              tabPanel("Read In Data", 
                                       wellPanel(
                                         fluidRow(
                                           column(4,
                                                  textInput(inputId = "project_name",
                                                            label = "Enter project name (no space) below:",
                                                            value = "Enter project name"))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_count",
                                                            label = "Select Count Table File",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_info",
                                                            label = "Select Design File",
                                                            accept = c(".csv"))))),
                                       
                                       
                                       wellPanel(
                                         tags$h4("View Count Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_count')),
                                       wellPanel(
                                         tags$h4("View Sample Information Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_info'))),
                              
                              
                              ## ------- exploratory analysis ---------
                              tabPanel("Exploratory Analysis", 
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(2,
                                                  selectInput(inputId = "gene_col_selected_expl",
                                                              label = "Select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(5,
                                                  selectInput(inputId = "sample_label_selected_expl",
                                                              label = "Select samples (same order as in histogram and boxplot):",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "covariate_selected_expl",
                                                              label = "Group by (select one):",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(2,
                                                  actionButton(inputId = "generate_result_expl",
                                                               label = "Generate Result")),
                                           column(2,
                                                  downloadButton(outputId = "download_expl",
                                                                 label = "Download Results"))
                                         )),
                                       wellPanel(
                                         tags$h3("Count Histogram"),
                                         plotlyOutput(outputId = "display_count_hist")),
                                       
                                       wellPanel(
                                         tags$h3("Between-sample Distribution: Boxplot"),
                                         fluidRow(
                                           column(6,
                                                  plotlyOutput(outputId = "display_boxplot")))),
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  wellPanel(
                                                    tags$h3("Heatmap Sampel-to-sample Distance"),
                                                    fluidRow(
                                                      plotOutput(outputId = "display_heatmap_expl")
                                                    ))),
                                           column(6,
                                                  wellPanel(
                                                    tags$h3("PCA Plot"),
                                                    fluidRow(plotlyOutput(outputId = "display_pca_expl"))))))),
                              
                              
                              ## ------------- DEGseq --------------
                              tabPanel("DE Analysis: DEGseq",
                                       # Data preparation
                                       wellPanel(tags$h3("Data Preparation"),
                                                 fluidRow(
                                                   column(2,
                                                          selectInput(inputId = "gene_col_selected_degseq",
                                                                      label = "Select gene column",
                                                                      choices = NULL,
                                                                      multiple = FALSE,
                                                                      selected = NULL)),
                                                   column(4,
                                                          selectInput(inputId = "sample_label_selected_degseq",
                                                                      label = "Select 2 Samples in order: trt1(denominator) trt2(numerator)",
                                                                      choices = NULL,
                                                                      multiple = TRUE,
                                                                      selected = NULL)),
                                                   column(4,
                                                          textInput(inputId = "row_sum",
                                                                    label = "Remove if sum across 2 samples is <",
                                                                    value = "10"))),
                                                 fluidRow(
                                                   column(2,
                                                          actionButton(inputId = "trim",
                                                                       label = "Generate Count Table"))),
                                                 fluidRow(
                                                   column(6,
                                                          tags$h4("Trimmed Count Table Summary:"),
                                                          verbatimTextOutput("trimmed_left"),
                                                          verbatimTextOutput("display_trimmed_ct")))),
                                       
                                       ## DE analysis
                                       wellPanel(
                                         tags$h3("DE Analysis"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "degseq_q_value",
                                                              label = "q-value(Storey et al. 2003) <",
                                                              choices = c(0.01, 0.05, 0.1, 0.25),
                                                              selected = 0.01)),
                                           column(3,
                                                  selectInput(inputId = "degseq_fold_change",
                                                              label = "Obs(log2(fold change)) >=",
                                                              choices = c(0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                              selected = 1)),
                                           column(3,
                                                  actionButton(inputId = "run_DEGexp",
                                                               label = "RUN DEGexp")),
                                           column(3,
                                                  downloadButton(outputId = "download_degseq",
                                                                 label = "Download Results"))),
                                         wellPanel(tags$h4("Result Table:"),
                                           DT::dataTableOutput(outputId = "result_table1_degseq")),
                                         
                                         wellPanel(tags$h4("MA Plot:"),
                                                   fluidRow(column(3,
                                                                   textInput(inputId = "ma1_title",
                                                                             label = "Enter MA plot title below:",
                                                                             value = ""))),
                                                   fluidRow(column(6,
                                                                   plotlyOutput(outputId = "ma1_degseq")),
                                                            column(6,
                                                                   tableOutput(outputId = "degseq_sign_number")))))),
                              
                              
                              ## --------------- DEseq2 ----------------
                              tabPanel("DE Analysis: DEseq2",
                                       
                                       # Data Preparation
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(2,
                                                  selectInput(inputId = "gene_col_selected_deseq2",
                                                              label = "Select gene column",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  selectInput(inputId = "sample_label_selected_deseq2",
                                                              label = "Select samples for DE analysis",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(3,
                                                  textInput(inputId = "row_sum_deseq2",
                                                            label = "Remove if total across samples is <",
                                                            value = "10")),
                                           column(2,
                                                  selectInput(inputId = "covariate_selected_deseq2",
                                                              label = "Select covariate(s)",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(3,
                                                  radioButtons(inputId = "design_deseq2",
                                                               label = "Select Design Formula:",
                                                               choices = list("Simple Comparison" = "o_interaction",
                                                                              "Combine Level" = "combine_level"
                                                                              # "With Interaction" = "w_interaction", 
                                                               ),
                                                               selected = NULL)),
                                           column(2,
                                                  actionButton(inputId = "generate_dds_from_matrix",
                                                               label = "Generate DESeqDataSet From Matrix"))),
                                         fluidRow(
                                           column(6,
                                                  tags$h4("DESeqDataSet Summary:"),
                                                  verbatimTextOutput(outputId = "trim_left_number"),
                                                  verbatimTextOutput("display_formula")))),
                                       
                                       # DE analysis
                                       wellPanel(
                                         tags$h3("DE Analysis"),
                                         fluidRow(column(4,
                                                         selectInput(inputId = "deseq2_p_value",
                                                                     label = "FDR adjusted p-value <",
                                                                     choices = c(0.01, 0.05, 0.1, 0.25),
                                                                     selected = 0.01)),
                                                  column(4,
                                                         selectInput(inputId = "deseq2_fold_change",
                                                                     label = "Obs(log2(fold change)) >=",
                                                                     choices = c(0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                                     selected = 1)),
                                                  column(3,
                                                         actionButton(inputId = "run_deseq",
                                                                      label = "Run DESeq"))),
                                         fluidRow(column(4,
                                                         selectInput(inputId = "comp1",
                                                                     label = "Select numerator level for the fold change",
                                                                     choices = NULL,
                                                                     multiple = FALSE,
                                                                     selected = NULL)),
                                                  column(4,
                                                         selectInput(inputId = "comp2",
                                                                     label = "Select denominator leve for the fold changel",
                                                                     choices = NULL,
                                                                     multiple = FALSE,
                                                                     selected = NULL)),
                                                  column(2,
                                                         actionButton(inputId = "extract_results_deseq2",
                                                                      label = "Extract Results")),
                                                  column(2,
                                                         downloadButton(outputId = "download_deseq2",
                                                                        label = "Download Results"))),
                                         fluidRow(verbatimTextOutput(outputId = "deseq2_results_dir")),
                                         
                                         wellPanel(tags$h4("Result Table:"),
                                                   DT::dataTableOutput(outputId = "display_dtf_res_contrast")),
                                         wellPanel(tags$h4("MA Plot:"),
                                                   fluidRow(column(3,
                                                                   textInput(inputId = "enter_deseq2_ma_title",
                                                                             label = "Enter MA plot title below:",
                                                                             value = ""))),
                                                   fluidRow(column(6,
                                                                   plotOutput(outputId = "display_deseq2_ma")),
                                                            column(6,
                                                                   tableOutput(outputId = "deseq2_sign_number")))),
                                         wellPanel(tags$h4("Plot Counts"),
                                                   fluidRow(column(3,
                                                                   textInput(inputId = "enter_gene_name",
                                                                             label = "Enter one gene name:")),
                                                            column(3,
                                                                   selectInput(inputId = "intgroup",
                                                                               label = "Select interested covariate(s) :",
                                                                               choices = NULL,
                                                                               multiple = TRUE,
                                                                               selected = NULL)),
                                                            column(2,
                                                                   actionButton(inputId = "generate_count_plot",
                                                                                label = "Generate Count Plot"))),
                                                   fluidRow(column(6,
                                                                   plotOutput(outputId = "deseq2_count_plot")))))),
                              ## ----------- change in gene expression analysis----------------
                              tabPanel("Change in Gene Expression",
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_contrast_1",
                                                            label = "Select RNA File of contrast 1",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_contrast_2",
                                                            label = "Select RNA File of contrast 1",
                                                            accept = c(".csv"))))),
                                       wellPanel(
                                         tags$h4("View Contrast 1 Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_contrast_1')),
                                       wellPanel(
                                         tags$h4("View Contrast 2 Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_contrast_2')),
                                       
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "gene_col_selected_rna_contrast_1",
                                                              label = "From Contrast 1 select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "logfold_change_selected_contrast_1",
                                                              label = "Select (normalized) logfold change column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "p_value_selected_rna_contrast_1",
                                                              label = "Select P value column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "gene_col_selected_rna_contrast_2",
                                                              label = "From Contrast 2 select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "logfold_change_selected_contrast_2",
                                                              label = "Select (normalized) logfold change column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "p_value_selected_rna_contrast_2",
                                                              label = "Select P value column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "p_thresh_change_in_gene_expr",
                                                              label = "p value <",
                                                              choices = c(0.01, 0.05, 0.1, 0.25),
                                                              selected = 0.01)),
                                           column(3,
                                                  selectInput(inputId = "logfold_change_thresh_change_in_gene_expr",
                                                              label = "Obs(log2(fold change)) >=",
                                                              choices = c(0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                              selected = 1)),
                                           column(3,
                                                  actionButton(inputId = "inner_join_2contrasts_rna",
                                                               label = "Inner Join Contrast 1 and Contrast 2 table"))),
                                         fluidRow(
                                           column(7,
                                                  tags$h4("Inner Joined Table Summary:"),
                                                  verbatimTextOutput("display_summary_inner_joined_2contrasts_rna")))),
                                       wellPanel(
                                         tags$h3("Change-in-gene-expression Heatmap"),
                                         fluidRow(
                                           column(3,
                                                  textInput(inputId = "change_in_gene_expr_plot_title",
                                                            label = "Enter plot title",
                                                            value = "Change in Gene Expression")),
                                           
                                           column(2,
                                                  textInput(inputId = "change_in_gene_expr_plot_x_text_1",
                                                            label = "x axis text 1",
                                                            value = "trt_2 vs trt_1")),
                                           column(2,
                                                  textInput(inputId = "change_in_gene_expr_plot_x_text_2",
                                                            label = "x axis text 2",
                                                            value = "trt_3 vs trt_2")),
                                           column(2,
                                                  actionButton(inputId = "generate_change_in_gene_expr_plot",
                                                               label = "Generate Heatmap")),
                                           column(2,
                                                  downloadButton(outputId = "download_change_in_gene_expr_plot",
                                                                 label = "Download Results"))),
                                         hr(),
                                         fluidRow(
                                           column(5,
                                                  plotOutput(outputId = "display_venn_diagram1_rna")),
                                           column(5,
                                                  plotOutput(outputId = "display_venn_diagram2_rna"))),
                                         hr(),
                                         fluidRow(
                                           column(7,
                                                  plotlyOutput(outputId = "change_in_gene_expr_plot")))))
                              # tabItem "rna-seq_analysis" ends
                              )),
                    
                    ## ---------------- DNA methylseq analysis ----------------
                    tabItem(tabName = "dna-seq_analysis",
                            tabsetPanel(
                              ## ------------ annotation -----------
                              tabPanel("Annotation",
                                       wellPanel(tags$h3("Annotation"),
                                                 shinyFilesButton(id = "dir_anno",
                                                                  label = "Select Peak file",
                                                                  title = "Please select a file",
                                                                  multiple = FALSE,
                                                                  class = c("csv")),
                                                 tags$h4("File Directory"),
                                                 verbatimTextOutput(outputId = "display_dir_anno"),
                                                 fluidRow(column(2,
                                                                 textInput(inputId = "tssreg_from",
                                                                           label = "TSS range from",
                                                                           value = -3000)),
                                                          column(2,
                                                                 textInput(inputId = "tssreg_to",
                                                                           label = "TSS range to",
                                                                           value = 3000)),
                                                          column(2,
                                                                 textInput(inputId = "txdb",
                                                                           label = "TxDb object",
                                                                           value = "TxDb.Mmusculus.UCSC.mm10.knownGene")),
                                                          column(2,
                                                                 textInput(inputId = "annodb",
                                                                           label = "Annotation package",
                                                                           value = "org.Mm.eg.db"))),
                                                 actionButton(inputId = "annotate",
                                                              label = "Annotate"))),
                              ## ----- exploratory analysis -----------
                              tabPanel("Exploratory Analysis",
                                       wellPanel(fluidRow(column(3,
                                                                 selectInput(inputId = "dna_sample_cols1",
                                                                             label = "Select sample columns for dtN",
                                                                             choices = NULL,
                                                                             multiple = TRUE,
                                                                             selected = NULL)),
                                                          column(4,
                                                                 selectInput(inputId = "dna_sample_cols2",
                                                                             label = "Select sample columns for dtX (same order as dtN)",
                                                                             choices = NULL,
                                                                             multiple = TRUE,
                                                                             selected = NULL)),
                                                          
                                                          column(3,
                                                                 textInput(inputId = "sample_name",
                                                                           label = "Enter sample name in order, seperated by comma",
                                                                           value = "")),
                                                          column(2,
                                                                 actionButton(inputId = "generate_expl_dna",
                                                                              label = "Generate Results"))),
                                                 verbatimTextOutput("test_dna")),
                                       
                                       wellPanel(tags$h3("Annotation by Region (%)"),
                                                 plotOutput(outputId = "anno_by_reg")),
                                       wellPanel(tags$h3("CpG Histogram by Region"),
                                                 plotOutput(outputId = "cpg_hist")),
                                       wellPanel(tags$h3("Methyl% by region"),
                                                 plotOutput(outputId = "meth_by_region")),
                                       # heatmap
                                       wellPanel(tags$h3("Heatmap of Methyl% on Gene Level"),
                                                 fluidRow(
                                                   column(3,
                                                          selectInput(inputId = "dna_heat_N",
                                                          label = "Select sample columns for dtN",
                                                          choices = NULL,
                                                          multiple = TRUE,
                                                          selected = NULL)),
                                                   column(3,
                                                          selectInput(inputId = "dna_heat_X",
                                                          label = "Select sample columns for dtX",
                                                          choices = NULL,
                                                          multiple = TRUE,
                                                          selected = NULL)),
                                                   column(3,
                                                          textInput(inputId = "dna_heat_sample_name",
                                                          label = "Enter sample name in order, seperated by comma",
                                                          value = ""))),
                                                 fluidRow(
                                                   column(3,
                                                          selectInput(inputId = "dna_heat_reg",
                                                          label = "Select one region of interest",
                                                          choices = c("Promoter", "5' UTR", "Body", "3' UTR", "Downstream"),
                                                          multiple = FALSE,
                                                          selected = "Promoter")),
                                                   column(4,
                                                          textInput(inputId = "dna_heat_top_n",
                                                          label = "Display top n genes with largest positive and negative mythly% diff",
                                                          value = "20")),
                                                   column(3,
                                                          actionButton(inputId = "generate_dna_heat",
                                                          label = "Generate Result"))),
                                                 fluidRow(
                                                   column(6,
                                                          plotlyOutput(outputId = "dna_heat1")),
                                                   column(6,
                                                          plotlyOutput(outputId = "dna_heat2"))))),
                              
                              ## ------------ DSS ---------------
                              tabPanel("DE Analysis: DSS",
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         tags$h4("Instruction"),
                                         p("To make BSseq Data for pairwise comparison, select the sample columns for each comparison, following the order 
                                           (1)N: Read coverage of the position from BS-seq data; 
                                           (2)X: Number of reads showing methylation of the position.",
                                           style = "padding-left: 0em"),
                                         p("For example, comparison1 (Con_N, Con_X), comparison2 (Exp_N, Exp_X). 
                                           Corresponding name can be Con and  Exp. For more information, please refer to  ",
                                           a("DSS. ", href = "https://bioconductor.org/packages/release/bioc/manuals/DSS/man/DSS.pdf"),
                                           style = "padding-left: 0em"),
                                         
                                         
                                         tags$h4("Construct BSseqData"),
                                         fluidRow(
                                           
                                           column(3,
                                                  selectInput(inputId = "dss_comp1",
                                                              label = "Select sample columns for comparison1",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "dss_comp2",
                                                              label = "Select sample columns for comparison2",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(3,
                                                  textInput(inputId = "dss_comp1_name",
                                                            label = "Enter sample name for comparison1",
                                                            value = "")),
                                           column(3,
                                                  textInput(inputId = "dss_comp2_name",
                                                            label = "Enter sample name for comparison2",
                                                            value = ""))
                                         ),
                                         tags$h4("Set Parameters for DML test"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "dss_wo_rep",
                                                              label = "Sample without replicates",
                                                              choices = c(TRUE,FALSE),
                                                              multiple = FALSE,
                                                              selected = TRUE)),
                                           column(3,
                                                  selectInput(inputId = "dss_smoothing",
                                                              label = "Smoothing",
                                                              choices = c(TRUE,FALSE),
                                                              multiple = FALSE,
                                                              selected = TRUE)),
                                           column(3,
                                                  textInput(inputId = "dss_smoothing_span",
                                                            label = "smoothing.span",
                                                            value = 500))),
                                         tags$h4("Set threshold"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "dss_pvalue",
                                                              label = "P value <",
                                                              choices = c(0.01, 0.05, 0.1, 0.25),
                                                              selected = 0.01)),
                                           column(3,
                                                  selectInput(inputId = "dss_methyl_diff",
                                                              label = "Obs(diff methyl%) >=",
                                                              choices = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
                                                              selected = 0.1)),
                                           column(2,
                                                  actionButton(inputId = "generate_dss",
                                                               label = "Generate result and download")))),
                                         wellPanel(
                                         tags$h3("Result of DML test"),
                                         DT::dataTableOutput(outputId = "dss_dml_tb"))),
                              
                              ## ------------ change in methyl ratio analysis ---------------
                              tabPanel("Change in Methylation Ratio",
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "dna_contrast_1",
                                                            label = "Select DNA File of contrast 1",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "dna_contrast_2",
                                                            label = "Select DNA File of contrast 1",
                                                            accept = c(".csv"))))),
                                       wellPanel(
                                         tags$h4("View Contrast 1 Table:"),
                                         DT::dataTableOutput(outputId = 'display_contrast_1')),
                                       wellPanel(
                                         tags$h4("View Contrast 2 Table:"),
                                         DT::dataTableOutput(outputId = 'display_contrast_2')),
                                       
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "gene_col_selected_dna_contrast_1",
                                                              label = "From Contrast 1 table select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "methyl_diff_dna_contrast_1",
                                                              label = "Select DNA methylation ratio difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "gene_col_selected_dna_contrast_2",
                                                              label = "From Contrast 2 table select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "methyl_diff_dna_contrast_2",
                                                              label = "Select DNA methylation ratio difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  actionButton(inputId = "inner_join_2contrasts",
                                                               label = "Inner Join Contrast 1 and Contrast 2 table"))),
                                         fluidRow(
                                           column(7,
                                                  tags$h4("Inner Joined Table Summary:"),
                                                  verbatimTextOutput("display_summary_inner_joined_2contrasts")))),
                                       wellPanel(
                                         tags$h3("Change-in-Methyl-Ratio Heatmap"),
                                         fluidRow(
                                           column(3,
                                                  textInput(inputId = "change_in_methyl_plot_title",
                                                            label = "Enter plot title",
                                                            value = "Change in Methylation Ratio")),
                                           
                                           column(2,
                                                  textInput(inputId = "change_in_methyl_plot_x_text_1",
                                                            label = "x axis text 1",
                                                            value = "trt_2 vs trt_1")),
                                           column(2,
                                                  textInput(inputId = "change_in_methyl_plot_x_text_2",
                                                            label = "x axis text 2",
                                                            value = "trt_3 vs trt_2")),
                                           column(1,
                                                  actionButton(inputId = "generate_change_in_methyl_plot",
                                                               label = "Generate Plot")),
                                           column(1,
                                                  downloadButton(outputId = "download_change_in_methyl_plot",
                                                                 label = "Download Results"))),
                                         hr(),
                                         fluidRow(
                                           column(7,
                                                  plotlyOutput(outputId = "change_in_methyl_plot")))))
                              # tabItem DNA end
                              )),
                    
                    ## --------------- DNA vs RNA -----------------
                    tabItem(tabName = "dna_vs_rna",
                            tabsetPanel(
                              tabPanel("Read In data",
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_sig",
                                                            label = "Select RNA File of Significant Genes",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "dna_sig",
                                                            label = "Select DNA File of Significant Genes",
                                                            accept = c(".csv"))))),
                                       wellPanel(
                                         tags$h4("View RNA Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_sig')),
                                       wellPanel(
                                         tags$h4("View DNA Table:"),
                                         DT::dataTableOutput(outputId = 'display_dna_sig'))),
                              tabPanel("Plot Data",
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "rna_gene_col_selected_rna_vs_dna",
                                                              label = "From RNA table select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "logfoldchange_col_selected_rna_vs_dna",
                                                              label = "Select RNA expression difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "dna_gene_col_selected_rna_vs_dna",
                                                              label = "From DNA table select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "methyl_diff_col_selected_rna_vs_dna",
                                                              label = "Select DNA methylation ratio difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           
                                           column(3,
                                                  selectInput(inputId = "region_col_selected_rna_vs_dna",
                                                              label = "Select region column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(2,
                                                  actionButton(inputId = "inner_join_rna_vs_dna",
                                                               label = "Inner Join RNA and DNA table"))),
                                         fluidRow(
                                           column(3,
                                                  textInput(inputId = "region_name_contian_rna_vs_dna",
                                                            label = "(Optional) Rename region that contain:")),
                                           column(3,
                                                  textInput(inputId = "region_new_name_rna_vs_dna",
                                                            label = "As:")),
                                           column(2,
                                                  actionButton(inputId = "rename_region_rna_vs_dna",
                                                               label = "Apply"))),
                                         fluidRow(
                                           column(7,
                                                  tags$h4("Inner Joined Table Summary:"),
                                                  verbatimTextOutput("display_summary_inner_joined_rna_dna")))),
                                       wellPanel(
                                         tags$h3("Starburst Plot"),
                                         fluidRow(
                                           column(3,
                                                  numericInput(inputId = "starburst_rna_hline",
                                                               label = "Enter RNA expression diff threshold",
                                                               value = 0.1)),
                                           column(3,
                                                  numericInput(inputId = "starburst_dna_vline",
                                                               label = "Enter DNA methylation % diff threshold",
                                                               value = 10)),
                                           column(3,
                                                  textInput(inputId = "starburst_title",
                                                            label = "Enter plot title",
                                                            value = "DNA vs RNA: starburst plot")),
                                           column(1,
                                                  actionButton(inputId = "generate_starburst_plot_rna_vs_dna",
                                                               label = "Generate Plot")),
                                           column(1,
                                                  downloadButton(outputId = "download_starburst_plot_and_joined_table",
                                                                 label = "Download Results"))),
                                         hr(),
                                         fluidRow(
                                           column(7,
                                                  plotlyOutput(outputId = "starburst_plot")))))
                            # tabItem DNA vs RNA end
                            ))
                    )))

## ------------- sssssSERVERrrrrr -------------------
server <- function(input, output, session) {
  ## ------------- upstream ---------------------
  # Specify folder containing FastQ file.
  # NOTE: all the files in the folder will be passed to the pipeline
  rv <- reactiveValues()
  
  shinyDirChoose(input = input, 
                 id = "directory", 
                 roots = volumes, 
                 session = session, 
                 restrictions = system.file(package = "base"))
  
  
  output$fileDir <- renderPrint({
    parseDirPath(roots = volumes, 
                 selection = input$directory)
  })
  
  output$program <- renderPrint({
    rv$pgmrna <- paste0("sbatch --job-name=",
                        input$job_name,
                        " --ntasks=",
                        input$ntasks,
                        " --cpus-per-taks=",
                        input$cpus_per_task,
                        " --mem=",
                        input$mem,
                        " --time=",
                        input$time,
                        " --output=",
                        input$output,
                        " --error=",
                        input$error,
                        ifelse(input$partition == "main", 
                               "",
                               paste0(" --partition=",
                                      input$partition)),
                        " module use /projects/community/modulefiles", 
                        " module load nextflow",
                        " export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)",
                        " export NXF_OPTS=\'-Xms1g -Xmx4g\'", 
                        " export NF_Work_Dir=\"/scratch/${USER}/NFWorkDir/${PWD}/work\"",
                        " export IS_SINGLE=\"",
                        input$isSingleEnd,
                        "\"",
                        " export IN_FILES=\"",
                        ifelse(input$isSingleEnd == "true",
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*.gz\""),
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*_R{1,2}*.fastq.gz\"")),
                        " mkdir -p $NF_Work_Dir",
                        " date",
                        " SECONDS=0",
                        " srun nextflow run rna-seq.nf --SingleEnd=$IS_SINGLE --reads=$IN_FILES -w $NF_Work_Dir -with-trace -with-report ${SLURM_JOB_PARTITION}_${SLURM_JOB_NODELIST}_${SLURM_JOB_ID}_RNA-seq-nf.html \\-with-timeline ${SLURM_JOB_PARTITION}_${SLURM_JOB_NODELIST}_${SLURM_JOB_ID}_RNA-seq-nf-timeline.html -resume",
                        " date",
                        " duration=$SECONDS",
                        " echo \"$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed.\"",
                        " touch \"${SLURM_JOB_PARTITION}-${SLURM_JOB_NODELIST}-${SLURM_JOB_ID}-$(($duration / 3600 ))h_$(($(($duration % 3600)) / 60 ))m_$(($duration % 60))s.time\""
    )
    
    rv$pgmdna <- paste0("sbatch --job-name=",
                        input$job_name,
                        " --ntasks=",
                        input$ntasks,
                        " --cpus-per-taks=",
                        input$cpus_per_task,
                        " --mem=",
                        input$mem,
                        " --time=",
                        input$time,
                        " --output=",
                        input$output,
                        " --error=",
                        input$error,
                        ifelse(input$partition == "main", 
                               "",
                               paste0(" --partition=",
                                      input$partition)),
                        " module use /projects/community/modulefiles", 
                        " module load nextflow",
                        " export NXF_OPTS=\'-Xms1g -Xmx4g\'",
                        " export NF_Work_Dir=\"/scratch/${USER}/NFWorkDir/${PWD}/work\"",
                        " mkdir -p $NF_Work_Dir",
                        " IS_SINGLE=\"",
                        input$isSingleEnd,
                        "\"",
                        " IN_FILES=\"",
                        ifelse(input$isSingleEnd == "true",
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*.gz\""),
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*_R{1,2}*.fastq.gz\"")),
                        # add the code inbetween?
                        " mkdir -p $NF_Work_Dir",
                        " date",
                        " SECONDS=0",
                        " srun nextflow run methyl-seq.nf -w $NF_Work_Dir --SingleEnd=$IS_SINGLE --reads=$IN_FILES -with-trace -with-report DNA-methyl.html  -with-timeline DNA-methyl-timeline.html -resume",
                        " EXIT_STATUS=$?",
                        " date",
                        " echo \"$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed.\"",
                        " touch \"${SLURM_JOB_PARTITION}-${SLURM_JOB_NODELIST}-${SLURM_JOB_ID}-$(($duration / 3600 ))h_$(($(($duration % 3600)) / 60 ))m_$(($duration % 60))s.time\"",
                        " if  [ $EXIT_STATUS  -eq  0 ]; then echo \" Pipeline completed. Removing WorkDir files\" rm -rf $NF_Work_Dir else echo \"Pipeline not completed.\" fi" 
    )
    
    cat(
      ifelse(input$is_rna == TRUE,
             paste('RNA Program:',
                   rv$pgmrna,
                   sep = '\n'),
             paste('DNA Program:',
                   rv$pgmdna, 
                   sep = '\n'))
    )
  })
  
  observeEvent(input$align, {
    # trigger nextflow to run rv$pgmrna or rv$pgmdna
    
  })
  
  ## ------------  RNA seq: Read in data and update SelectInput -------------
  
  observeEvent(input$rna_count,{
    infile_rna_count <- input$rna_count
    tb <- fread(infile_rna_count$datapath, 
                skip = 1)
    rv$ct <- as.data.frame(tb)
    
    output$display_rna_count = DT::renderDataTable({
      DT::datatable(rv$ct,
                    options = list(scrollX = TRUE))})
  })
  
  observeEvent(input$rna_info,{
    infile_rna_info <- input$rna_info
    tb <- fread(infile_rna_info$datapath, 
                header = TRUE)
    tb <- as.data.frame(tb)
    rv$info <- tb
    
    output$display_rna_info = DT::renderDataTable({
      DT::datatable(rv$info,
                    options = list(scrollX = TRUE))})
    
    updateSelectInput(session, 
                      inputId = "gene_col_selected_expl", 
                      label = "Select gene column", 
                      choices = colnames(rv$ct), 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "sample_label_selected_expl", 
                      label = "Select samples (same order as in histogram and boxplot)", 
                      choices = rv$info$`Sample Label`, 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "covariate_selected_expl",
                      label = "Group by (select one):",
                      choices = colnames(rv$info)[c(-1,-2)], 
                      selected = NULL)
    
    updateSelectInput(session, 
                      inputId = "gene_col_selected_degseq", 
                      label = "Select gene column", 
                      choices = colnames(rv$ct), 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "sample_label_selected_degseq", 
                      label = "Select 3 samples in order: trt1 trt2 trt3", 
                      choices = rv$info$`Sample Label`, 
                      selected = NULL)
    
    updateSelectInput(session, 
                      inputId = "gene_col_selected_deseq2", 
                      label = "Select gene column", 
                      choices = colnames(rv$ct), 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "sample_label_selected_deseq2", 
                      label = "Select samples for DE analysis", 
                      choices = rv$info$`Sample Label`, 
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "covariate_selected_deseq2",
                      label = "Select Covariate(s)",
                      choices = colnames(rv$info)[c(-1,-2)], 
                      selected = NULL)
    
  })
  
  ## ---------------------- Exploratory Analysis -------------
  observeEvent(input$generate_result_expl,{
    
    withProgress(message = "Generating Result",
                 value = 0,
                 expr = {
                   # prepare selected count table
                   incProgress(0.2, detail = "Preparing Data")
                   rv$ct[is.na(rv$ct)] <- 0
                   dt <- try(
                     match.ct.info(rv$ct, 
                                   rv$info, 
                                   input$sample_label_selected_expl, 
                                   input$gene_col_selected_expl), 
                     silent = T)
                   if (class(dt) == "try-error") {
                     showNotification(dt[1],
                                      type = "error",
                                      duration = 15)
                     
                   }else{
                     # selected design table
                     selected_info <- rv$info[match(input$sample_label_selected_expl, 
                                                    rv$info$`Sample Label`), 
                                              c("Sample Label", 
                                                input$covariate_selected_expl)]
                     colnames(selected_info) <- c("sample", "covariate")
                     # log transform
                     dt_log2 <- dt
                     dt_log2[, 2:ncol(dt_log2)] <- log2(dt_log2[, 2:ncol(dt_log2)] + 1)
                     dt_log2_long <- melt(dt_log2,
                                          id.vars = 1,
                                          measure.vars = 2:ncol(dt_log2),
                                          variable.name = "sample",
                                          value.name = "count")
                     # merged long 
                     dt1 <- merge(dt_log2_long, selected_info, by = "sample")
                     
                     # histogram
                     rv$count_hist <- ggplot(data = dt1,
                                             aes(x = count)) +
                       geom_histogram() +
                       xlab("") +
                       ylab("log2(count + 1)") + 
                       facet_wrap( ~ sample)
                     
                     p1 <- ggplotly(rv$count_hist)
                     output$display_count_hist <- renderPlotly({print(p1)})
                     
                     # boxplot
                     rv$boxplot <- ggplot(data = dt1,
                                          aes(x = sample,
                                              y = count,
                                              fill = covariate )) + 
                       geom_boxplot() +
                       xlab("") +
                       ylab("log2(count + 1)") + 
                       labs(fill = "legend")
                     
                     p2 <- ggplotly(rv$boxplot)
                     output$display_boxplot <- renderPlotly({print(p2)})
                     
                     # heatmap s2s
                     tmp <- dt_log2[, -1]
                     sample_dist <- as.matrix(dist(t(as.matrix(tmp))))

                     rv$expl_heatmap <- pheatmap(sample_dist,
                                                 col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))

                     output$display_heatmap_expl <- renderPlot({print(rv$expl_heatmap)})

                     # pca
                     matrix_pca <- t(dt_log2[, -1])
                     matrix_pca <- matrix_pca[, apply(matrix_pca, 2, var) != 0]
                     res_pca <- prcomp(matrix_pca,
                                       center = TRUE,
                                       scale. = TRUE)
                     df_out <- as.data.frame(res_pca$x)

                     percentage <- round(res_pca$sdev^2 / sum(res_pca$sdev^2) * 100, 2)
                     percentage <- paste( colnames(df_out), "(", paste0( as.character(percentage), "%"), ")")

                     df_out$sample <- colnames(dt_log2)[-1]
                     df_out <- merge(df_out, selected_info, by = "sample")

                     rv$pca <- ggplot(data = df_out,
                                      aes(x = PC1,
                                          y = PC2)) +
                       geom_point(aes(fill = df_out[, ncol(df_out)]),
                                  shape = 21,
                                  size = 3,
                                  alpha = 0.5) +
                       xlab(percentage[1]) +
                       ylab(percentage[2]) +
                       geom_text(data = df_out,
                                 aes(x = PC1 + 20,
                                     y = PC2,
                                     label = df_out$sample),
                                 size = 2,
                                 hjust = 0.5) +
                       labs(fill = "legend")

                     p4 <- ggplotly(rv$pca)

                     output$display_pca_expl <- renderPlotly({print(p4)})
                     
                     # end of if else   
                   }
                   
                   
                   incProgress(0.8, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
    
  })
  
  ## ---------------------- DEGseq  ----------------
  ## Data Preparation
  observeEvent(input$trim,{
    if (length(input$sample_label_selected_degseq) != 2) {
      showNotification("Please select two samples in order",
                       type = "error",
                       duration = 15)
    }else{
      dt <- try(
        match.ct.info(rv$ct, 
                      rv$info, 
                      input$sample_label_selected_degseq, 
                      input$gene_col_selected_degseq), 
        silent = T)
      if (class(dt) == "try-error") {
        showNotification(dt[1],
                         type = "error",
                         duration = 15)
        
      }else{
        # trim 
        rv$ct_degseq <- subset(dt,
                               rowSums(dt[, -1]) > as.numeric(input$row_sum))
        
        output$trimmed_left <- renderPrint({
          cat(
            paste(nrow(rv$ct_degseq),
                  "genes left, down from",
                  nrow(rv$ct),
                  "genes",
                  sep = " "))})
        output$display_trimmed_ct <- renderPrint({summary(rv$ct_degseq)})
      }                         
    }  
  })
  
  ## DE analysis
  observeEvent(input$run_DEGexp,{
    withProgress(message = "Running DEGexp",
                 value = 0,
                 expr = {
                   # run DEGexp
                   incProgress(0.4, 
                               detail = "Processing trt2-trt1")
                   outputDir <-  file.path(tempdir())
                   rv$ct_degseq_matrix <- as.matrix(rv$ct_degseq)
                   DEGexp(geneExpMatrix1 = rv$ct_degseq_matrix,
                          geneCol1 = 1, 
                          expCol1 = 3, 
                          groupLabel1 = colnames(rv$ct_degseq_matrix)[3],
                          
                          geneExpMatrix2 = rv$ct_degseq_matrix,
                          geneCol2 = 1, 
                          expCol2 = 2,
                          groupLabel2 = colnames(rv$ct_degseq_matrix)[2],
                          
                          foldChange = as.numeric(input$degseq_fold_change),
                          
                          qValue = as.numeric(input$degseq_q_value),
                          thresholdKind = 5, 
                          rawCount = TRUE,
                          normalMethod = "none",
                          method = "MARS",
                          outputDir = outputDir)
                   
                   rv$table_trt2_trt1 <- fread(paste(outputDir,"\\output_score.txt", sep = ""))
                   
                   # add mu clr and pch
                   incProgress(0.2, detail = "Displaying Result Tables")
                   
                   rv$table_trt2_trt1 <- add.clm(rv$table_trt2_trt1, 
                                                 input$degseq_q_value, 
                                                 input$degseq_fold_change)
                   rv$table_trt2_trt1_round <- rv$table_trt2_trt1
                   rv$table_trt2_trt1_round[, c(4:9,11)] <- round(rv$table_trt2_trt1_round[, c(4:9,11)], 2)
                   
                   
                   
                   
                   output$result_table1_degseq = DT::renderDataTable({
                     DT::datatable(rv$table_trt2_trt1_round,
                                   filter = "top",
                                   options = list(scrollX = TRUE,
                                                  columnDefs = list(list(searchable = FALSE, 
                                                                         targets = c(1:4,6:7,11:13)))))
                   })
                   
                   # MA plot
                   incProgress(0.2, detail = "Generating MA plot")
                   
                   rv$ma1_degseq <- ma(rv$table_trt2_trt1, input$degseq_q_value,
                                       input$degseq_fold_change) 
                   p1 <- ggplotly(rv$ma1_degseq)
                   output$ma1_degseq <- renderPlotly({print(p1)})
                   
                   result <- c("up-regulated DEG",
                                 "non-significant DEG",
                                 "down-regulated DEG")
                     number <- c(sum(rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$degseq_q_value) & 
                                       rv$table_trt2_trt1$`log2(Fold_change) normalized` >= as.numeric(input$degseq_fold_change), 
                                     na.rm = TRUE),
                                 nrow(rv$table_trt2_trt1) - sum(rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$degseq_q_value) & 
                                                                   abs(rv$table_trt2_trt1$`log2(Fold_change) normalized`) >= as.numeric(input$degseq_fold_change),
                                                                 na.rm = TRUE),
                                 sum(rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$degseq_q_value) & 
                                       rv$table_trt2_trt1$`log2(Fold_change) normalized` <= -as.numeric(input$degseq_fold_change),
                                     na.rm = TRUE))
                     degseq_sign_number_tb <- data.frame(result, number)
                     output$degseq_sign_number <- renderTable(degseq_sign_number_tb)
                   
                   
                   
                   
                   incProgress(0.2, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
  })
  
  # renew MA plot title
  observeEvent(input$ma1_title,{
    if (!is.null(rv$ma1_degseq)) {
      rv$ma1_degseq <- rv$ma1_degseq +
        ggtitle(input$ma1_title)
      p1 <- ggplotly(rv$ma1_degseq)
      
      output$ma1_degseq <- renderPlotly({print(p1)})
    }
  })
  
  # observeEvent(input$ma2_title,{
  #   if (!is.null(rv$ma2_degseq)) {
  #     rv$ma2_degseq <- rv$ma2_degseq +
  #       ggtitle(input$ma2_title) 
  #     p2 <- ggplotly(rv$ma2_degseq)
  #     output$ma2_degseq <- renderPlotly({print(p2)})
  #   }
  # })
  
  
  
  ## ---------------------- DEseq2 -----------------------
  observeEvent(input$generate_dds_from_matrix,{
    rv$sample_name_list_deseq2 <- rv$info$`Sample Name`[match(input$sample_label_selected_deseq2, 
                                                              rv$info$`Sample Label`)]
    rv$ct_deseq2 <- rv$ct[, rv$sample_name_list_deseq2]
    rv$ct_deseq2 <- as.matrix(rv$ct_deseq2)
    rownames(rv$ct_deseq2) <- rv$ct[, input$gene_col_selected_deseq2]
    colnames(rv$ct_deseq2) <- input$sample_label_selected_deseq2
    rv$ct_deseq2[is.na(rv$ct_deseq2)] <- 0
    
    
    rv$info_deseq2 <- rv$info[match(input$sample_label_selected_deseq2,
                                    rv$info$`Sample Label`), 
                              c("Sample Label",
                                input$covariate_selected_deseq2)]
    colnames(rv$info_deseq2) <- c("sample", 
                                  input$covariate_selected_deseq2)
    rv$info_deseq2[-1] <- lapply(rv$info_deseq2[-1], factor)
    
    updateSelectInput(session, 
                      inputId = "intgroup", 
                      label = "Select interested covariate(s):", 
                      choices = colnames(rv$info_deseq2)[-1], 
                      selected = NULL)
    
    if (length(input$covariate_selected_deseq2) == 1 &
       input$design_deseq2 == "o_interaction") {
      rv$formula <- as.formula(paste("~",
                                     input$covariate_selected_deseq2))
    }else if (length(input$covariate_selected_deseq2) > 1 &
             input$design_deseq2 == "combine_level") {
      rv$info_deseq2$group <- factor(apply(rv$info_deseq2[ , input$covariate_selected_deseq2],
                                           1 ,
                                           paste ,
                                           collapse = "_" ))
      rv$formula <- as.formula("~ group")
    }else{
      rv$formula <- "incorrect"
    } 
    
    output$display_formula <- renderPrint({
      cat(
        paste0("The design formula: ", 
               paste(as.character(rv$formula), 
                     collapse = " "))
      )
    })
    
    if (rv$formula == "incorrect") {
      showNotification("Incorrect design formula",
                       type = "error",
                       duration = 15)
    }else{
      rv$dds <- try(DESeqDataSetFromMatrix(countData = rv$ct_deseq2,
                                           colData = rv$info_deseq2,
                                           design = rv$formula),
                    silent = TRUE)
      if (class(rv$dds) == "try-error") {
        showNotification(rv$dds[1],
                         type = "error",
                         duration = 15)
        output$trim_left_number <- renderPrint({
          cat(rv$dds)
        })
        
      }else{
        rv$dds_trimmed <- rv$dds[rowSums(counts(rv$dds)) >= as.numeric(input$row_sum_deseq2), ]
        output$trim_left_number <- renderPrint({
          cat(
            paste(nrow(counts(rv$dds_trimmed)),
                  "genes left, down from",
                  nrow(counts(rv$dds)),
                  "genes"))
        })
      } 
    }
  })
  
  observeEvent(input$run_deseq,{
    
    if (!is.null(rv$dds_trimmed)) {
      snowparam <- SnowParam(workers = snowWorkers(), 
                             type = "SOCK")
      register(snowparam, 
               default = TRUE)
      
      withProgress(message = "Running DESeq", 
                   value = 0,
                   expr = {
                     incProgress(0.2, detail = "Processing 20%")
                     incProgress(0.3, detail = "Processing 50%")
                     incProgress(0.2, detail = "Processing 70%")
                     rv$dds_res <- DESeq(rv$dds_trimmed,
                                         fitType = "local",
                                         parallel = FALSE) # TRUE
                     incProgress(0.3, detail = "Processing 100%")
                     Sys.sleep(1)
                   })
      
      # update contrast choice
      if (input$design_deseq2 == "o_interaction") {
        contrast <- levels(rv$info_deseq2[,2])
      }else if (input$design_deseq2 == "combine_level") {
        contrast <- levels(rv$info_deseq2$group)
      }
      updateSelectInput(session,
                        inputId = "comp1",
                        label = "Select numerator level for the fold change",
                        choices = contrast,
                        selected = NULL)
      
      updateSelectInput(session,
                        inputId = "comp2",
                        label = "Select denominator level for the fold change",
                        choices = contrast,
                        selected = NULL)
      
    }else{
      showNotification("please generate dds object from matrix first",
                       type = "error",
                       duration = 15)
    }  
  })
  
  observeEvent(input$extract_results_deseq2,{
    if (!is.null(rv$dds_res)) {
      withProgress(message = "Extracting results with contrast", 
                   value = 0,
                   expr = {
                     # extract result
                     incProgress(0.2, detail = "Preparing result table")
                     if (input$design_deseq2 == "o_interaction") {
                       contrast <- c(input$covariate_selected_deseq2, input$comp1, input$comp2)
                     }else if (input$design_deseq2 == "combine_level") {
                       contrast <- c("group", input$comp1, input$comp2)
                     }
                     rv$dds_res_contrast <- results(rv$dds_res,
                                                    contrast = contrast)
                     rv$dtf_res_contrast <- data.frame(gene = rownames(rv$dds_res_contrast),
                                                       do.call("cbind", 
                                                               rv$dds_res_contrast@listData))
                     
                     rv$dtf_res_contrast$clr <- "black"
                     rv$dtf_res_contrast$clr[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value)] <- "purple"
                     rv$dtf_res_contrast$clr[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & 
                                               abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change)] <- "red"
                     rv$dtf_res_contrast$pch <- 46
                     rv$dtf_res_contrast$pch[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value)] <- 3
                     rv$dtf_res_contrast$pch[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & 
                                               abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change)] <- 4
                     rv$dtf_res_contrast$pch <- factor(rv$dtf_res_contrast$pch,
                                                       levels = c(46, 3, 4))
                     
                     rv$dtf_res_contrast_round <- rv$dtf_res_contrast
                     rv$dtf_res_contrast_round[, c(2:7)] <- round(rv$dtf_res_contrast_round[, c(2:7)], 2)
                     
                     output$display_dtf_res_contrast <- DT::renderDataTable({
                       DT::datatable(rv$dtf_res_contrast_round,
                                     filter = "top",
                                     options = list(scrollX = TRUE,
                                                    columnDefs = list(list(searchable = FALSE, 
                                                                           targets = c(1:2, 4:6, 8:9)))))
                     })
                     
                     # MA plot
                     incProgress(0.2, detail = "Preparing MA plot")
                     rv$deseq2_ma <- ggplot(rv$dtf_res_contrast,
                                            aes(x = log(baseMean + 1),
                                                y = `log2FoldChange`,
                                                colour = clr,
                                                shape = pch)) +
                       scale_shape_manual(name = "Legend:",
                                          values = c(46, 3, 4),
                                          labels = c("No significance",
                                                     paste("p-Value < ", 
                                                           as.numeric(input$deseq2_p_value)),
                                                     paste("p-Value < ", 
                                                           as.numeric(input$deseq2_p_value), 
                                                           " & abs(log2) >= ", 
                                                           as.numeric(input$deseq2_fold_change)))) +
                       scale_color_manual(name = "Legend:",
                                          values = c("black",
                                                     "purple",
                                                     "red"),
                                          labels = c("No significance",
                                                     paste("p-Value < ", 
                                                           as.numeric(input$deseq2_p_value)),
                                                     paste("p-Value < ", 
                                                           as.numeric(input$deseq2_p_value), 
                                                           " & abs(log2) >= ", 
                                                           as.numeric(input$deseq2_fold_change)))) +
                       geom_hline(yintercept = c(-as.numeric(input$deseq2_fold_change), 
                                                 as.numeric(input$deseq2_fold_change)),
                                  lty = 2) +
                       theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank(),
                             axis.line = element_line(colour = "black"),
                             plot.title = element_text(hjust = 0.5),
                             legend.position = "top") + 
                       geom_point() +
                       labs(x = "log2(BaseMean + 1)", y = "log2FoldChange")
                     output$display_deseq2_ma <- renderPlot({
                       print(rv$deseq2_ma)
                     })
                     
                     # sign gene number
                     result <- c("up-regulated DEG",
                                 "non-significant DEG",
                                 "down-regulated DEG")
                     number <- c(sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & 
                                       rv$dtf_res_contrast$log2FoldChange >= as.numeric(input$deseq2_fold_change), 
                                     na.rm = TRUE),
                                 nrow(rv$dtf_res_contrast) - sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & 
                                                                   abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change),
                                                                 na.rm = TRUE),
                                 sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & 
                                       rv$dtf_res_contrast$log2FoldChange <= -as.numeric(input$deseq2_fold_change),
                                     na.rm = TRUE))
                     deseq2_sign_number_tb <- data.frame(result, number)
                     output$deseq2_sign_number <- renderTable(deseq2_sign_number_tb)
                     
                     incProgress(0.6, detail = "Processing 100%")
                     Sys.sleep(2)
                   })
      # renew MA plot
      observeEvent(input$enter_deseq2_ma_title,{
        rv$deseq2_ma <- rv$deseq2_ma + 
          ggtitle(input$enter_deseq2_ma_title)
        # p1 <- ggplotly(rv$deseq2_ma)
        output$display_deseq2_ma <- renderPlot({
          print(rv$deseq2_ma)
        })
      })
    }
  })
  
  ## generate count plot
  observeEvent(input$generate_count_plot,{
    if (is.null(input$enter_gene_name) | is.null(input$intgroup)) {
      showNotification("Please enter one gene name and select interested covariate(s)",
                       type = "error",
                       duration = 15)
    }else if (input$enter_gene_name %in% rownames(rv$dds_trimmed)  == FALSE ) {
      showNotification("Please enter the exact correct gene name",
                       type = "error",
                       duration = 15)
    }else{
      dt <- plotCounts(rv$dds_res, gene = input$enter_gene_name,
                       intgroup = input$intgroup, 
                       returnData = TRUE)
      if (length(input$intgroup) > 1) {
        dt$group <- apply(dt[, -1], 
                          1, 
                          paste,
                          collapse = "_")
      }
      
      rv$count_plot <- ggplot(dt,
                              aes(x = dt[,ncol(dt)],
                                  y = count)) +
        geom_point(position = position_jitter(width = 0.1, height = 0)) +
        labs(y = "normalized count", 
             x = NULL)
      output$deseq2_count_plot  <- renderPlot({print(rv$count_plot)})
    }
  })
  
  ## ---------------------- RNA: download results ----------------
  output$download_expl <- downloadHandler(
    filename = function() {
      paste("exploratory_results", 
            ".zip", 
            sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_expl_count_hist",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 6,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$count_hist)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_expl_boxplot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path2,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$boxplot)
      dev.off()
      
      file_path3 <- paste(input$project_name,
                          "_expl_pca",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path3,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$pca)
      dev.off()
      
      file_path4 <- paste(input$project_name,
                          "_expl_heatmap",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path4,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$expl_heatmap)
      dev.off()
      
      fs <- c(file_path1, file_path2, file_path3, file_path4)
      zip(zipfile = fname, files = fs)
    },
    contentType = "application/zip"
  )
  
  
  output$download_degseq <- downloadHandler(
    filename = function() {
      paste("degseq_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_degseq_ma1",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$ma1_degseq)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "degseq_trt2_trt1.csv",
                          sep = "")
      write.csv(rv$table_trt2_trt1, file_path2, row.names = FALSE)
      
      fs <- c(file_path1,
              file_path2)
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
  

  
  output$download_deseq2 <- downloadHandler(
    filename = function() {
      paste("deseq2_results", 
            ".zip", 
            sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_deseq2_result_table",
                          ".csv",
                          sep = "")
      write.csv(rv$dtf_res_contrast,
                file_path1,
                row.names = FALSE)
      
      
      file_path2 <- paste(input$project_name,
                          "_deseq2_ma",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path2,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$deseq2_ma)
      dev.off()
      
      if (!is.null(rv$count_plot)) {
        file_path4 <- paste(input$project_name,
                            "_deseq2_count_plot.tiff",
                            sep = "")
        tiff(filename = file_path4,
             height = 6,
             width = 6,
             units = 'in',
             res = 300,
             compression = "lzw+p")
        print(rv$count_plot)
        dev.off()
      }else{
        file_path4 <- NULL
      }
      
      fs <- c(file_path1, 
              file_path2,
              file_path4)
      zip(zipfile = fname,
          files = fs)
    },
    contentType = "application/zip"
  )
  
  ## ----------- change in gene expre analysis ------------
  observeEvent(input$rna_contrast_1,{
    infile_rna_contrast_1 <- input$rna_contrast_1
    tb <- read.table(infile_rna_contrast_1$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_rna_contrast_1 <- as_tibble(tb)
    output$display_rna_contrast_1 = DT::renderDataTable({
      DT::datatable(rv$tb_rna_contrast_1,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_rna_contrast_1",
                      label = "From Contrast 1 select gene column:",
                      choices = colnames(rv$tb_rna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "logfold_change_selected_contrast_1",
                      label = "Select (normalized) logfold change column:",
                      choices = colnames(rv$tb_rna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "p_value_selected_rna_contrast_1",
                      label = "Select P value column:",
                      choices = colnames(rv$tb_rna_contrast_1),
                      selected = NULL)})
  
  observeEvent(input$rna_contrast_2,{
    infile_rna_contrast_2 <- input$rna_contrast_2
    tb <- read.table(infile_rna_contrast_2$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_rna_contrast_2 <- as_tibble(tb)
    output$display_rna_contrast_2 = DT::renderDataTable({
      DT::datatable(rv$tb_rna_contrast_2,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_rna_contrast_2",
                      label = "From Contrast 2 select gene column:",
                      choices = colnames(rv$tb_rna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "logfold_change_selected_contrast_2",
                      label = "Select (normalized) logfold change column:",
                      choices = colnames(rv$tb_rna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "p_value_selected_rna_contrast_2",
                      label = "Select P value column:",
                      choices = colnames(rv$tb_rna_contrast_2),
                      selected = NULL)})
  
  observeEvent(input$inner_join_2contrasts_rna,{
    contrast1 <- rv$tb_rna_contrast_1 %>%
      filter( !!as.name(input$p_value_selected_rna_contrast_1) < input$p_thresh_change_in_gene_expr &
              abs(!!as.name(input$logfold_change_selected_contrast_1)) >= input$logfold_change_thresh_change_in_gene_expr  ) %>% 
      select(input$gene_col_selected_rna_contrast_1,
             input$logfold_change_selected_contrast_1) %>%
      rename("gene" = input$gene_col_selected_rna_contrast_1,
             "log_foldchange_1" = input$logfold_change_selected_contrast_1)
      
    
    contrast2 <- rv$tb_rna_contrast_2 %>%
      filter( !!as.name(input$p_value_selected_rna_contrast_2) < input$p_thresh_change_in_gene_expr &
              abs(!!as.name(input$logfold_change_selected_contrast_2)) >= input$logfold_change_thresh_change_in_gene_expr  ) %>% 
      select(input$gene_col_selected_rna_contrast_2,
             input$logfold_change_selected_contrast_2) %>%
      rename("gene" = input$gene_col_selected_rna_contrast_2,
             "log_foldchange_2" = input$logfold_change_selected_contrast_2)
    
    rv$joined_2contrasts_rna <- try(
      inner_join(contrast1, contrast2, by = "gene") %>%
        filter((log_foldchange_1*log_foldchange_2) < 0),
      silent = T)
    
    if (inherits(rv$joined_2contrasts_rna, "try-error")) {
      showNotification("Please choose the correct columns and try again",
                       type = "error",
                       duration = 15)

    }else{
      rv$rna_contrast1_up_num <-  nrow(filter(contrast1, log_foldchange_1 > 0))
      rv$rna_contrast1_down_num <-  nrow(filter(contrast1, log_foldchange_1 < 0))
      rv$rna_contrast2_up_num <-  nrow(filter(contrast2, log_foldchange_2 > 0))
      rv$rna_contrast2_down_num <-  nrow(filter(contrast2, log_foldchange_2 < 0))
      rv$rna_up_donw <- nrow(filter(rv$joined_2contrasts_rna, log_foldchange_1 > 0))
      rv$rna_down_up <- nrow(filter(rv$joined_2contrasts_rna, log_foldchange_1 < 0))
      output$display_summary_inner_joined_2contrasts_rna <- renderPrint({summary(rv$joined_2contrasts_rna)})
    }
    
    })
  
  observeEvent(input$generate_change_in_gene_expr_plot,{
    output$display_venn_diagram1_rna <- renderPlot({
      p1 <- draw.pairwise.venn(area1 = rv$rna_contrast1_up_num,
                               area2 = rv$rna_contrast2_down_num,
                               cross.area = rv$rna_up_donw,
                               scaled = TRUE,
                               col = c("green3", "firebrick"))
      rv$venn_diagram1_rna <- grid.arrange(gTree(children = p1), top = "Trt2-Trt1 Up Trt3-Trt2 Down")
    })
      
    output$display_venn_diagram2_rna <- renderPlot({
      p2 <- draw.pairwise.venn(area1 = rv$rna_contrast1_down_num,
                                       area2 = rv$rna_contrast2_up_num,
                                       cross.area = rv$rna_down_up,
                                       scaled = TRUE,
                                       col = c("firebrick", "green3"))
      rv$venn_diagram2_rna <- grid.arrange(gTree(children = p2), top = "Trt2-Trt1 Down Trt3-Trt2 Up")
    })
    rv$change_in_gene_expr_heatmap <-  two_column_heatmap(rv$joined_2contrasts_rna,
                                                          gene,
                                                          log_foldchange_1,
                                                          log_foldchange_2,
                                                          input$change_in_gene_expr_plot_title,
                                                          input$change_in_gene_expr_plot_x_text_1,
                                                          input$change_in_gene_expr_plot_x_text_2)
    
    p <- ggplotly(rv$change_in_gene_expr_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_gene_expr_plot <- renderPlotly({print(p)})
  })
  
  
  
  output$download_change_in_gene_expr_plot <- downloadHandler(
    filename = function() {
      paste("change_in_gene_expr_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_change_in_gene_expr_plot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 10,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$change_in_gene_expr_heatmap)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_change_in_gene_expr_table",
                          ".csv",
                          sep = "")
      write.csv(rv$joined_2contrasts_rna, 
                file_path2, 
                row.names = FALSE)
      file_path3 <- paste(input$project_name,
                          "_venn_up_donw.tiff",
                          sep = "")
      tiff(filename = file_path3,
           height = 4,
           width = 5,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      grid.draw(rv$venn_diagram1_rna)
      dev.off()
      
      file_path4 <- paste(input$project_name,
                          "_venn_donw_up.tiff",
                          sep = "")
      tiff(filename = file_path4,
           height = 4,
           width = 5,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      grid.draw(rv$venn_diagram2_rna)
      dev.off()
      
      fs <- c(file_path1, 
              file_path2,
              file_path3,
              file_path4)
      
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
  
  ## ------------  DNA mettyl: Annotate and update selectinput------------
  shinyFileChoose(input = input, 
                  id = "dir_anno", 
                  roots = volumes,
                  session = session,
                  filetypes = c('csv'))
  
  output$display_dir_anno <- renderPrint({
    parseFilePaths(
      roots = volumes, 
      selection = input$dir_anno)$datapath
  })
  
  
  observeEvent(input$annotate,
               withProgress(message = "Running annotatepeak",
                            value = 0,
                            expr = {
                              incProgress(0.7, detail = "Processing 70%")
                              rv$peakAnno1 <- annotatePeak(peak = parseFilePaths(roots = volumes,selection = input$dir_anno)$datapath, 
                                                           tssRegion = c(as.numeric(input$tssreg_from), 
                                                                         as.numeric(input$tssreg_to)), 
                                                           TxDb = eval(parse(text = input$txdb)),
                                                           annoDb = input$annodb)
                              
                              rv$dt_dna <- data.frame(start = rv$peakAnno1@anno@ranges@start,
                                                      as.data.frame(rv$peakAnno1@anno@elementMetadata@listData))
                              
                              
                              updateSelectInput(session, 
                                                inputId = "dna_sample_cols1", 
                                                label = "Select sample columns for dtN", 
                                                choices = colnames(rv$dt_dna), 
                                                selected = NULL)
                              updateSelectInput(session, 
                                                inputId = "dna_sample_cols2", 
                                                label = "Select sample columns for dtX (same order as dtN)", 
                                                choices = colnames(rv$dt_dna), 
                                                selected = NULL)
                              
                              updateSelectInput(session, 
                                                inputId = "dna_heat_N", 
                                                label = "Select sample columns for dtN", 
                                                choices = colnames(rv$dt_dna), 
                                                selected = NULL)
                              updateSelectInput(session, 
                                                inputId = "dna_heat_X", 
                                                label = "Select sample columns for dtX (same order as dtN)", 
                                                choices = colnames(rv$dt_dna), 
                                                selected = NULL)
                              updateSelectInput(session, 
                                                inputId = "dss_comp1", 
                                                label = "Select sample columns for comparison1", 
                                                choices = colnames(rv$dt_dna), 
                                                selected = NULL)
                              updateSelectInput(session, 
                                                inputId = "dss_comp2", 
                                                label = "Select sample columns for comparison2", 
                                                choices = colnames(rv$dt_dna), 
                                                selected = NULL)
                              
                              incProgress(0.3, detail = "Processing 100%")
                            })
  )
  
  ## ----------------------- Explotary analysis ------------------------
  observeEvent(input$generate_expl_dna,{
    
    # anno by reg (not related to samples)
    t1 <- rv$peakAnno1@annoStat
    t1$Feature <- factor(t1$Feature,
                         levels = as.character(t1$Feature[order(t1$Frequency,
                                                                decreasing = TRUE)]))
    
    rv$anno_by_reg <- ggplot(t1,
                             aes(x = "",
                                 y = Frequency,
                                 fill = Feature)) +
      geom_bar(width = 1, 
               stat = "identity",
               color = "black") +
      coord_polar("y",
                  start = 0,
                  direction = -1) +
      scale_x_discrete("") +
      ggtitle("Annotation by Region (%)")
    
    output$anno_by_reg <- renderPlot({
      print(rv$anno_by_reg)
    })
    
    # construct dt with samples selected by users
    # remove genes with unmapped region
    rv$dt_dna1 <- rv$dt_dna[!is.na(rv$dt_dna$SYMBOL == "NA"), ]
    
    rv$dt_dna2 <- data.frame(gene = rv$dt_dna1$SYMBOL,
                             anno = rv$dt_dna1$annotation,
                             geneId = rv$dt_dna1$geneId,
                             chr = rv$dt_dna1$geneChr,
                             pos = rv$dt_dna1$start,
                             reg = NA,
                             CpG = rv$dt_dna1$CpG,
                             rv$dt_dna1[, input$dna_sample_cols1],
                             rv$dt_dna1[, input$dna_sample_cols2],
                             geneName = rv$dt_dna1$GENENAME)
    
    
    rv$dt_dna2 <- try(droplevels(rv$dt_dna2[rowSums(rv$dt_dna2[, c(input$dna_sample_cols1, input$dna_sample_cols2)],
                                                    na.rm = TRUE) > 0, ]),
                      silent = TRUE)
    
    if (class(rv$dt_dna2) == "try-error") {
      showNotification(rv$dds[1],
                       type = "error",
                       duration = 15)
    }else{
      output$test_dna <- renderPrint( cat(paste0("Removed ",
                                                 nrow(rv$dt_dna) - nrow(rv$dt_dna1),
                                                 " genes with unmapped region and ",
                                                 nrow(rv$dt_dna1) - nrow(rv$dt_dna2),
                                                 " genes with all NA value, with ",
                                                 nrow(rv$dt_dna2), 
                                                 " genes left.")))
      rv$dt_dna2$reg <- as.character(rv$dt_dna2$anno)
      
      rv$dt_dna2$reg[substr(rv$dt_dna2$anno,
                            1,
                            8) == "Promoter"] <- "Promoter"
      rv$dt_dna2$reg[substr(rv$dt_dna2$anno,
                            1,
                            4) %in% c("Exon",
                                      "Intr")] <- "Body"
      rv$dt_dna2$reg[substr(rv$dt_dna2$anno,
                            1,
                            4) %in% c("Dist",
                                      "Down")] <- "Downstream"
      rv$dt_dna2$reg <- factor(rv$dt_dna2$reg,
                               levels = c("Promoter",
                                          "5' UTR",
                                          "Body",
                                          "3' UTR",
                                          "Downstream"))
      
      # cpg hist (not related to samples)
      rv$cpg_hist <- ggplot(rv$dt_dna2,
                            aes(x = CpG)) +
        facet_wrap(~ reg, scale = "free_y") +
        geom_histogram(color = "black",
                       fill = "grey",
                       binwidth = 5) +
        scale_x_continuous(name = "Number of CpG",
                           breaks = c(3,
                                      seq(from = 5,
                                          to = 60,
                                          by = 5))) +
        coord_cartesian(xlim = c(3, 60)) +
        scale_y_continuous(name = "Counts") +
        ggtitle("Distribution of DMR by Number of CpG and Region")
      
      output$cpg_hist <- renderPlot({
        print(rv$cpg_hist)
      })
      
      
      # avergae methy% by region by sample
      dtN <- as.matrix(rv$dt_dna2[, input$dna_sample_cols1])
      dtX <- as.matrix(rv$dt_dna2[, input$dna_sample_cols2])
      # Add 0.5 to all NAs and zeros in meth. hits
      # NOTE: if there were no hits (N = NA or 0), the pct will be NA anyway
      dtX <- apply(dtX,
                   2,
                   function(a) {
                     a[is.na(a)] <- a[a == 0] <- 0.5
                     return(a)
                   })
      pct <- dtX/dtN
      
      dt.pct <- data.table(rv$dt_dna2$reg,
                           pct)
      
      colnames(dt.pct) <- c("reg",
                            unlist(strsplit(input$sample_name, ",")))
      
      dt.pct.mean <- aggregate( . ~ reg,
                                dt.pct, mean)
      
      dt.pct.long <- melt(dt.pct.mean,
                          id.vars = "reg",
                          variable.name = "trt",
                          value.name = "Methylation (%)")
      
      output$test_dna1 <- renderPrint(dt.pct.long)
      output$test_dna2 <- renderPrint(dt.pct.mean)
      rv$meth_by_region <- ggplot(dt.pct.long,
                                  aes(x = reg,
                                      y = 100*`Methylation (%)`,
                                      group = trt,
                                      fill = trt)) +
        geom_bar(position = position_dodge(),
                 stat = "identity") +
        scale_x_discrete("Region") +
        scale_y_continuous("Average Methylation (%)",
                           limits = c(0, 100)) +
        scale_fill_discrete("Treatment") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45,
                                         hjust = 1))
      
      output$meth_by_region <- renderPlot({
        print(rv$meth_by_region)
      })  
    }
  })
  
  # Heatmap of methyl% on gene level
  observeEvent(input$generate_dna_heat, {
    dt <- rv$dt_dna[!is.na(rv$dt_dna$SYMBOL == "NA"), ]
    
    dt1 <- data.frame(gene = dt$SYMBOL,
                      anno = dt$annotation,
                      reg = NA,
                      dt[, input$dna_heat_X],
                      dt[, input$dna_heat_N],
                      geneName = dt$GENENAME)
    
    dt1$reg <- as.character(dt1$anno)
    rv$dt_dna2$reg <- as.character(rv$dt_dna2$anno)
    dt1$reg[substr(dt1$anno,
                   1,
                   8) == "Promoter"] <- "Promoter"
    dt1$reg[substr(dt1$anno,
                   1,
                   4) %in% c("Exon",
                             "Intr")] <- "Body"
    dt1$reg[substr(dt1$anno,
                   1,
                   4) %in% c("Dist",
                             "Down")] <- "Downstream"
    dt1$reg <- factor(dt1$reg,
                      levels = c("Promoter",
                                 "5' UTR",
                                 "Body",
                                 "3' UTR",
                                 "Downstream"))
    dt2 <- na.omit(dt1)
    
    dtN <- as.matrix(dt2[, input$dna_heat_N])
    dtX <- as.matrix(dt2[, input$dna_heat_X])
    pct <- as.data.frame(dtX/dtN)
    colnames(pct) <- unlist(strsplit(input$dna_heat_sample_name, ","))
    pct$gene <- dt2$gene
    pct$reg <- dt2$reg
    # filter reg
    pct <- pct[pct$reg == input$dna_heat_reg, ]
    pct$diff_trt2_trt1 <- pct[, 2] - pct[, 1] 
    pct_sort <- pct[order(pct$diff_trt2_trt1, 
                          decreasing = TRUE) , ]  # trt1,trt2,trt3,gene,diff
    pct_sort$diff_trt2_trt1 <- NULL # trt1,trt2,trt3,gene,reg
    pct_sort$reg <- NULL
    
    pct_sort_pos <- pct_sort[1:as.numeric(input$dna_heat_top_n) ,]
    pct_sort_neg <- pct_sort[(nrow(pct_sort) - 1 - as.numeric(input$dna_heat_top_n)):nrow(pct_sort), ]
    
    pos_long <- melt(data = pct_sort_pos,
                     id.vars = "gene",
                     measure.vars = unlist(strsplit(input$dna_heat_sample_name, ",")))
    neg_long <- melt(data = pct_sort_neg,
                     id.vars = "gene",
                     measure.vars = unlist(strsplit(input$dna_heat_sample_name, ",")))
    
    rv$dna_heat1 <- ggplot(data = pos_long) +
      geom_tile(aes(x = variable,
                    y = gene,
                    fill = value),
                color = "black") +
      scale_fill_gradient2(high = "red",
                           limit = c(0, 1),
                           name = "Methylation") +
      scale_x_discrete("Treatment",
                       expand = c(0, 0)) +
      scale_y_discrete("Gene",
                       expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 30,
                                       hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0("Top ", 
                     input$dna_heat_top_n,
                     " Genes With Largest\nPositive Differences in ",
                     colnames(pct_sort_pos)[2], 
                     " vs ",
                     colnames(pct_sort_pos)[1], 
                     " in ", 
                     input$dna_heat_reg ))
    
    output$dna_heat1 <- renderPlotly({
      print(ggplotly(rv$dna_heat1))
    })
    
    rv$dna_heat2 <- ggplot(data = neg_long) +
      geom_tile(aes(x =  variable,
                    y = gene,
                    fill = value),
                color = "black") +
      scale_fill_gradient2(high = "red",
                           limit = c(0, 1),
                           name = "Methylation") +
      scale_x_discrete("Treatment",
                       expand = c(0, 0)) +
      scale_y_discrete("Gene",
                       expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 30,
                                       hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0("Top ", 
                     input$dna_heat_top_n,
                     "  Genes With Largest\nNegative Differences in ",
                     colnames(pct_sort_pos)[2], 
                     " vs ",
                     colnames(pct_sort_pos)[1], 
                     " in ", 
                     input$dna_heat_reg ))
    
    output$dna_heat2 <- renderPlotly({
      print(ggplotly(rv$dna_heat2))
    })
    
    
  })
  
  ## ----------------------- DSS -------------
  fdtdss <- reactive({
    dt <- data.frame(start = rv$peakAnno1@anno@ranges@start,
                     as.data.frame(rv$peakAnno1@anno@elementMetadata@listData))
    dt <- dt[!is.na(dt$SYMBOL == "NA"), ]
    
    comp1 <- input$dss_comp1
    comp2 <- input$dss_comp2
    name1 <- input$dss_comp1_name
    name2 <- input$dss_comp2_name
    
    df1 <- data.frame(chr = dt$geneChr,
                      pos = dt$start,
                      N = dt[, comp1[1]],
                      X = dt[, comp1[2]])
    df2 <- data.frame(chr = dt$geneChr,
                      pos = dt$start,
                      N = dt[, comp2[1]],
                      X = dt[, comp2[2]])
    frames <- list(df1, df2)
    names(frames) <- c(name1, name2)
    bsdata <- makeBSseqData(frames, names(frames))
    
    # perform DML test, differntially methylated loci (DML) for two group comparisons of bisulfite sequencing (BS-seq)
    # without replicates, must set equal.disp=T
    
    dml <- DMLtest(bsdata, 
                   group1 = name1,
                   group2 = name2,
                   equal.disp = TRUE,
                   smoothing = TRUE, 
                   smoothing.span = 500)
    dml
  })
  
  observeEvent(input$generate_dss, {
    dml <- fdtdss()

    write.csv(dml,
              file = paste("dss_",
                           Sys.Date(),
                           ".csv",
                           sep = ""),
              row.names = FALSE)
    
    output$dss_dml_tb <- DT::renderDataTable({
      datatable(dml,
                    options = list(scrollX = TRUE))})
  })
  
  # # This should work but doesn't
  # # Source: https://yihui.shinyapps.io/DT-info/
  # output$download_dss <- downloadHandler(
  #   filename = function() {
  #     paste("dss_",
  #           Sys.Date(),
  #           ".csv",
  #           sep = "")
  #   },
  #   content = function(file) {
  #     dml <- input$dss_dml_tb_rows_all
  #     write.csv(dml,
  #               file = file,
  #               row.names = FALSE)
  #   })
  
  ## --------------- Change-in-Methyl Heatmap -------------
  observeEvent(input$dna_contrast_1,{
    infile_dna_contrast_1 <- input$dna_contrast_1
    tb <- read.table(infile_dna_contrast_1$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_dna_contrast_1 <- as_tibble(tb)
    output$display_contrast_1 = DT::renderDataTable({
      DT::datatable(rv$tb_dna_contrast_1,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_dna_contrast_1",
                      label = "From Contrast 1 table select gene column:",
                      choices = colnames(rv$tb_dna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "methyl_diff_dna_contrast_1",
                      label = "Select DNA methylation ratio difference column:",
                      choices = colnames(rv$tb_dna_contrast_1),
                      selected = NULL)})
  
  observeEvent(input$dna_contrast_2,{
    infile_dna_contrast_2 <- input$dna_contrast_2
    tb <- read.table(infile_dna_contrast_2$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_dna_contrast_2 <- as_tibble(tb)
    output$display_contrast_2 = DT::renderDataTable({
      DT::datatable(rv$tb_dna_contrast_2,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_dna_contrast_2",
                      label = "From Contrast 2 table select gene column:",
                      choices = colnames(rv$tb_dna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "methyl_diff_dna_contrast_2",
                      label = "Select DNA methylation ratio difference column:",
                      choices = colnames(rv$tb_dna_contrast_2),
                      selected = NULL)})
  
  observeEvent(input$inner_join_2contrasts,{
    contrast1 <- rv$tb_dna_contrast_1 %>%
      select(input$gene_col_selected_dna_contrast_1,
             input$methyl_diff_dna_contrast_1) %>%
      rename("gene" = input$gene_col_selected_dna_contrast_1,
             "methyl_diff_1" = input$methyl_diff_dna_contrast_1)
    
    contrast2 <- rv$tb_dna_contrast_2  %>%
      select(input$gene_col_selected_dna_contrast_2,
             input$methyl_diff_dna_contrast_2)  %>%
      rename("gene" = input$gene_col_selected_dna_contrast_2,
             "methyl_diff_2" = input$methyl_diff_dna_contrast_2)
    
    rv$joined_2contrasts <- try(
      inner_join(contrast1, contrast2, by = "gene") %>% 
        filter((methyl_diff_1*methyl_diff_2) < 0),
      silent = T)
    
    if (inherits(rv$joined_2contrasts, "try-error")) {
      showNotification("Please choose the correct columns and try again",
                       type = "error",
                       duration = 15)

    }else{
      output$display_summary_inner_joined_2contrasts <- renderPrint({summary(rv$joined_2contrasts)})
    }
    
    })
  
  observeEvent(input$generate_change_in_methyl_plot,{
    rv$change_in_methyl_heatmap <-  two_column_heatmap(rv$joined_2contrasts,
                                                             gene,
                                                             methyl_diff_1,
                                                             methyl_diff_2,
                                                             input$change_in_methyl_plot_title,
                                                             input$change_in_methyl_plot_x_text_1,
                                                             input$change_in_methyl_plot_x_text_2)
    
    p <- ggplotly(rv$change_in_methyl_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_methyl_plot <- renderPlotly({print(p)})
  })
  
  
  
  output$download_change_in_methyl_plot <- downloadHandler(
    filename = function() {
      paste("change_in_methyl_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_change_in_methyl_plot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 10,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$change_in_methyl_heatmap)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_change_in_methyl_table",
                          ".csv",
                          sep = "")
      write.csv(rv$joined_2contrasts, 
                file_path2, 
                row.names = FALSE)
      
      fs <- c(file_path1,
              file_path2)
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
  
  ## -------------  RNA vs DNA: Read in data and update SelectInput ------------
  
  observeEvent(input$rna_sig,{
    infile_rna_sig <- input$rna_sig
    tb <- read.table(infile_rna_sig$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_rna_sig <- as_tibble(tb)
    output$display_rna_sig = DT::renderDataTable({
      DT::datatable(rv$tb_rna_sig,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "rna_gene_col_selected_rna_vs_dna",
                      label = "From RNA table select gene column:",
                      choices = colnames(rv$tb_rna_sig),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "logfoldchange_col_selected_rna_vs_dna",
                      label = "Select RNA expression difference column:",
                      choices = colnames(rv$tb_rna_sig),
                      selected = NULL)})
  
  observeEvent(input$dna_sig,{
    infile_dna_sig <- input$dna_sig
    tb <- read.table(infile_dna_sig$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_dna_sig <- as_tibble(tb)
    output$display_dna_sig = DT::renderDataTable({
      DT::datatable(rv$tb_dna_sig,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "dna_gene_col_selected_rna_vs_dna",
                      label = "From DNA table select gene column:",
                      choices = colnames(rv$tb_dna_sig),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "methyl_diff_col_selected_rna_vs_dna",
                      label = "Select DNA methylation ratio difference column:",
                      choices = colnames(rv$tb_dna_sig),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "region_col_selected_rna_vs_dna",
                      label = "Select region column:",
                      choices = colnames(rv$tb_dna_sig),
                      selected = NULL)})
  
  observeEvent(input$inner_join_rna_vs_dna,{
    rna_sig <- rv$tb_rna_sig %>%
      select(input$rna_gene_col_selected_rna_vs_dna,
             input$logfoldchange_col_selected_rna_vs_dna) %>%
      rename("gene" = input$rna_gene_col_selected_rna_vs_dna,
             "RNA_exp_diff" = input$logfoldchange_col_selected_rna_vs_dna)
    
    dna_sig <- rv$tb_dna_sig  %>%
      select(input$dna_gene_col_selected_rna_vs_dna,
             input$methyl_diff_col_selected_rna_vs_dna,
             input$region_col_selected_rna_vs_dna)  %>%
      rename("gene" = input$dna_gene_col_selected_rna_vs_dna,
             "DNA_methyl_diff" = input$methyl_diff_col_selected_rna_vs_dna,
             "region" = input$region_col_selected_rna_vs_dna)
    
    rv$joined_rna_dna <- try(
      inner_join(rna_sig, dna_sig, by = "gene"),
      silent = T)
    
    if (inherits(rv$joined_rna_dna, "try-error")) {
      showNotification("Please choose the correct columns and try again",
                       type = "error",
                       duration = 15)

    }else{
      rv$joined_rna_dna %>%
        mutate(region = as.factor(region)) -> rv$joined_rna_dna
      output$display_summary_inner_joined_rna_dna <- renderPrint({summary(rv$joined_rna_dna)})
    }
    
    })
  
  
  observeEvent(input$rename_region_rna_vs_dna,{
    rv$joined_rna_dna %>%
      mutate(region = as.character(region)) %>%
      mutate(region = replace(region,
                              str_detect(region,
                                         input$region_name_contian_rna_vs_dna),
                              input$region_new_name_rna_vs_dna)) %>%
      mutate(region =  as.factor(region)) -> 
      rv$joined_rna_dna
    output$display_summary_inner_joined_rna_dna <- renderPrint({summary(rv$joined_rna_dna)})

  })
  
  observeEvent(input$generate_starburst_plot_rna_vs_dna,{
    rv$starburst_plot <- starburst(rv$joined_rna_dna,
                                   gene,
                                   RNA_exp_diff,
                                   DNA_methyl_diff,
                                   region,
                                   input$starburst_rna_hline,
                                   input$starburst_dna_vline,
                                   input$starburst_title)
    p1 <- ggplotly(rv$starburst_plot,
                   tooltip = c("text", "x", "y")) %>%
      layout(height = 800, width = 800)
    output$starburst_plot <- renderPlotly({print(p1)})
  })
  
  
  output$download_starburst_plot_and_joined_table <- downloadHandler(
    filename = function() {
      paste("rna_vs_dna_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_rna_vs_dna_starburst_plot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 10,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$starburst_plot)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_rna_vs_dna_joined_table",
                          ".csv",
                          sep = "")
      write.csv(rv$joined_rna_dna, 
                file_path2, 
                row.names = FALSE)
      
      fs <- c(file_path1,
              file_path2)
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)