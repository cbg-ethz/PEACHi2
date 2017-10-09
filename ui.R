library(shiny)

shinyUI(navbarPage(title="PEACHi2", inverse=TRUE,
                   tabPanel("Home",
                            titlePanel("Patterns of Expression and Analysis of Clusters of HIV/Host interactions v2.0"), 
                            p("This resource allows querying time series transcriptome, proteome and phosphoproteome responses of SupT1 T cells to HIV infection."),
                            tags$br(),
                            sidebarLayout(
                              sidebarPanel( width = 3,
                                            strong("Select proteo-transcriptomic behavior:"),
                                            tags$br(),
                                            tags$br(),
                                            radioButtons("rna_choice", "RNA:",
                                                         c("Up", "Down", "Not DE", "All"), inline = TRUE, selected = "Up"),
                                            radioButtons("prot_choice", "Protein:",
                                                         c("Up", "Down", "Not DE", "All"), inline = TRUE, selected = "Up"),
                                            radioButtons("phospho_choice", "Phosphoprotein",
                                                         c("Up", "Down", "Not DE", "All"), inline = TRUE, selected = "All"),
                                            actionButton("select1", "Submit"),
                                            tags$hr(style="background-color: black; height: 1px; border: 0;"),
                                            strong("Insert a list of genes:"),
                                            tags$br(),
                                            tags$br(),
                                            tags$textarea(id = 'list_genes', placeholder = 'Insert list of genes here', width="100%", rows = 10, cols = 30, ""),
                                            #submitButton("Submit"),
                                            tags$br(),
                                            tags$br(),
                                            actionButton("select2", "Submit"),
                                            tags$hr(style="background-color: black; height: 1px; border: 0;"),
                                            strong("Download displayed data table:"),
                                            tags$br(),
                                            tags$br(),
                                            downloadButton('downloadData', 'Download data table')
                              ),
                              mainPanel(
                                plotOutput(outputId = "timeSeries"),
                                tabsetPanel(id = "data_tabs",
                                            tabPanel('RNA',
                                                     tabsetPanel( id = "tabs_rna",
                                                       tabPanel("Normalized HIV/Mock ratios", DT::dataTableOutput("rna_tab"), value = 1),
                                                       tabPanel("Normalized read counts (M=Mock, H=HIV)", DT::dataTableOutput("rna_count_tab"), value = 2)
                                                     )
                                            ),
                                            tabPanel('Protein',
                                                     tabsetPanel( id = "tabs_prot",
                                                       tabPanel("Normalized HIV/Mock ratios", DT::dataTableOutput("protein_tab"), value = 1),
                                                       tabPanel("Peptide counts", DT::dataTableOutput("prot_count_tab"), value = 2)
                                                     )
                                            ),
                                            tabPanel('Phosphoprotein',
                                                     DT::dataTableOutput("phospho_tab"))
                                )
                              )
                            )
                   ),
                   tabPanel("Help",
                        img(src="PEACHi_help.jpg", style="float:center", height = 235.75089*2.5, width = 741.52722*2.5 ) 
                   ),
                   tabPanel("About"
                   )
))



#img(src="peach.jpg", style="float:left", height = 40, width = 68 ), 