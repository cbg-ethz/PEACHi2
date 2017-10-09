library(shiny)
library(DT)

source("scripts/load_data.R")
source("scripts/plot_time_series.R")
source("scripts/retrieve_data.R")

# pre-defined initial configuration of the website display
initial_data = retrieve_genes("Up", "Up", "All", rna, prot, phospho)
initial_table = retrieve_tables(as.data.frame(initial_data), list(rna_r = rna$Up, prot_r = prot$Up, phospho_r = phospho$All), rna$counts, prot$counts)

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  data_to_plot = reactiveValues(data = initial_data)
  data_from_table = reactiveValues(rna_tab = initial_table$rna_tab, rna_count_tab = initial_table$rna_count_tab, 
                                   prot_tab = initial_table$prot_tab, prot_count_tab = initial_table$prot_count_tab, 
                                   phospho_tab = initial_table$phospho_tab)
  plot_tables = reactiveValues(rna_r = rna$Up, prot_r = prot$Up, phospho_r = phospho$All)
  proxy_rna_tab = dataTableProxy('rna_tab')
  proxy_prot_tab = dataTableProxy('rna_tab')
  proxy_phospho_tab = dataTableProxy('rna_tab')
  
  observeEvent(input$select2, {
    if(length(grep(",",input$list_genes))>0) {
      data_to_plot$data = lapply(strsplit(input$list_genes,split = ','), function(x) trim(x))
    } else {
      data_to_plot$data = lapply(strsplit(input$list_genes,split = '\n'), function(x) trim(x))
    }
    plot_tables$rna_r = rna$All; plot_tables$prot_r = prot$All; plot_tables$phospho_r = phospho$All
    r_table = retrieve_tables(as.data.frame(data_to_plot$data), plot_tables, rna$counts, prot$counts)
    data_from_table$rna_tab = r_table$rna_tab; data_from_table$rna_count_tab = r_table$rna_count_tab; data_from_table$prot_tab = r_table$prot_tab
    data_from_table$prot_count_tab = r_table$prot_count_tab; data_from_table$phospho_tab = r_table$phospho_tab
    updateTextInput(session, "list_genes", value = "")    
  })  
  
  observeEvent(input$select1, {
    # clear gene list and remove selected rows
    updateTextInput(session, "list_genes", value = "")
    proxy_rna_tab %>% selectRows(list()) #selectRows(NULL)
    proxy_prot_tab %>% selectRows(list()) #selectRows(NULL)
    proxy_phospho_tab %>% selectRows(list()) #selectRows(NULL)
    # pharse the radio buttons labels
    rna_choice = strsplit(input$rna_choice," ")[[1]][1]
    prot_choice = strsplit(input$prot_choice," ")[[1]][1]
    ph_choice = strsplit(input$phospho_choice," ")[[1]][1]
    
    # get the desired list of genes and corresponding tables
    if (input$rna_choice == "All" & input$prot_choice == "All" & input$phospho_choice == "All") {
      plot_tables$rna_r = rna$All; plot_tables$prot_r = prot$All; plot_tables$phospho_r = phospho$All
      data_to_plot$data = "all"
      r_table = retrieve_tables(data_to_plot$data, plot_tables, rna$counts, prot$counts)
      data_from_table$rna_tab = r_table$rna_tab; data_from_table$rna_count_tab = r_table$rna_count_tab; data_from_table$prot_tab = r_table$prot_tab
      data_from_table$prot_count_tab = r_table$prot_count_tab; data_from_table$phospho_tab = r_table$phospho_tab
    } else {
      data_to_plot$data = retrieve_genes(rna_choice, prot_choice, ph_choice, rna, prot, phospho)
      plot_tables$rna_r = rna[[rna_choice]]; plot_tables$prot_r = prot[[prot_choice]]; plot_tables$phospho_r = phospho[[ph_choice]]
      r_table = retrieve_tables(as.data.frame(data_to_plot$data), plot_tables, rna$counts, prot$counts)
      data_from_table$rna_tab = r_table$rna_tab; data_from_table$rna_count_tab = r_table$rna_count_tab; data_from_table$prot_tab = r_table$prot_tab
      data_from_table$prot_count_tab = r_table$prot_count_tab; data_from_table$phospho_tab = r_table$phospho_tab

    }
  }) 
  
  
  # Plot the genes
  output$timeSeries <- renderPlot({
    selected_rows = list(rna_selected = data_from_table$rna_tab[input$rna_tab_rows_selected,], 
                         prot_selected = data_from_table$prot_tab[input$protein_tab_rows_selected,],
                         ph_selected = data_from_table$phospho_tab[input$phospho_tab_rows_selected,])
    plot_trio_given_list(plot_tables$rna_r, plot_tables$prot_r, plot_tables$phospho_r,  rna$cutoff, prot$cutoff, phospho$cutoff,
                         as.data.frame(data_to_plot$data), c("#ffa366","#9bb8e4","#b8e600"), selected_rows)
    
  })
  
  
  # Display data tables
  output$rna_tab = DT::renderDataTable({
    if (is.null(data_to_plot$data)) return()
    datatable(data_from_table$rna_tab, selection = 'single', options = list(pageLength = 25))
  })
  
  output$rna_count_tab = DT::renderDataTable({
    if (is.null(data_to_plot$data)) return()
    datatable(data_from_table$rna_count_tab, selection = 'none', options = list(pageLength = 25))
  })
  
  output$protein_tab = DT::renderDataTable({ 
    if (is.null(data_to_plot$data)) return()
    datatable(data_from_table$prot_tab, selection = 'single', options = list(pageLength = 25))
  })
  
  output$prot_count_tab = DT::renderDataTable({
    if (is.null(data_to_plot$data)) return()
    datatable(data_from_table$prot_count_tab, selection = 'none', options = list(pageLength = 25))
  })
  
  output$phospho_tab = DT::renderDataTable({
    if (is.null(data_to_plot$data)) return()
    datatable(data_from_table$phospho_tab, selection = 'single', options = list(pageLength = 25))
  })
  
  
  # Download data tables
  output$downloadData <- downloadHandler(
    filename = function() {paste("Data_",input$data_tabs,".csv",sep="")},
    content = function(file) {
      if (input$data_tabs == "RNA" && input$tabs_rna == 1) {
        filename = function() { "Data_RNA_FC.csv" }
        write.csv(data_from_table$rna_tab, file)
      } else if (input$data_tabs == "RNA" && input$tabs_rna == 2) {
          write.csv(data_from_table$rna_count_tab, file)
      } else if (input$data_tabs == "Protein" && input$tabs_prot == 1) {
        write.csv(data_from_table$prot_tab, file)
      } else if (input$data_tabs == "Protein" && input$tabs_prot == 2) {
          write.csv(data_from_table$prot_count_tab, file)
      } else if (input$data_tabs == "Phosphoprotein") {
        write.csv(data_from_table$phospho_tab, file)
      }
    }
  )
  
  
  
  
  
})