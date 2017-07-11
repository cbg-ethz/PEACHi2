#################
# Script for retrieving the data table for the queried list of genes
#################

retrieve_genes = function(rna_choice, prot_choice, ph_choice, rna, prot, phospho) {
  # extract first keyword ("Up", "Down", "All")
  rna_choice = strsplit(rna_choice," ")[[1]][1]
  prot_choice = strsplit(prot_choice," ")[[1]][1]
  ph_choice = strsplit(ph_choice," ")[[1]][1]
  
  list_genes = list()
  i = 0
  # retrieve the intersection of the desired behaviors
  if (rna_choice != "All") {
    i = i + 1
    list_genes[[i]] = rna[[rna_choice]][,1]
  }
  if (prot_choice != "All") {
    i = i + 1
    list_genes[[i]] = prot[[prot_choice]][,1]
  }
  if (ph_choice != "All") {
    i = i + 1
    list_genes[[i]] = phospho[[ph_choice]][,1]
  }
  
  genes = list_genes[[i]]
  while (i>0) {
    genes = intersect(genes,list_genes[[i]])
    i = i - 1
  }
  return(genes)
}

retrieve_table = function(gene_list, df_in, tab_names) {
  pos_data = pos_ann = annotated_df = NULL
  if (sum(gene_list == "all")==1) {
    annotated_df = merge(gene_annotations, df_in, by=1, all.y = TRUE)
  } else {
    pos_data = which(toupper(df_in[,1]) %in% toupper(gene_list[,1]))
    pos_ann = which(toupper(gene_annotations[,1]) %in% toupper(gene_list[,1]))
    if (!is.null(pos_ann) & !is.null(pos_data)) {
      annotated_df = merge(gene_annotations[pos_ann, ], df_in[pos_data, ], by=1, all.y = TRUE)
    }
  }
  if (!is.null(annotated_df)) {
    annotated_df = annotated_df[,-2] # remove the ENSG id
    annotated_df = cbind(annotated_df[,c(1,2)], apply(annotated_df[,3:ncol(annotated_df)], c(1,2), function(x) suppressWarnings(as.numeric(format(x, digits = 3)))))
    colnames(annotated_df) = tab_names
  }
 
  return(annotated_df)
}

retrieve_tables = function(gene_list, plot_tables, rna_count, prot_count) {
  table_data = NULL
  if (nrow(as.data.frame(gene_list))==0) return(table_data)
  
  table_data$rna_tab = retrieve_table(gene_list, plot_tables$rna_r, c("Gene HGNC symbol","Gene description", "6h", "12h", "18h", "24h", "30h"))
  table_data$rna_count_tab = retrieve_table(gene_list, rna_count, c("Gene HGNC symbol","Gene description", "M_6h", "M_12h", "M_18h", "M_24h_A", "M_24h_B", "M_30h", 
                                                                                 "H_6h", "H_12h", "H_18h", "H_24h_A", "H_24H_B", "H_30h"))
  table_data$prot_tab = retrieve_table(gene_list, plot_tables$prot_r, c("Gene HGNC symbol","Gene description", "2h", "6h", "12h", "18h", "24h"))
  table_data$prot_count_tab = retrieve_table(gene_list, prot_count, c("Gene HGNC symbol","Gene description", "2h", "6h", "12h", "18h", "24h"))
  table_data$phospho_tab = retrieve_table(gene_list, plot_tables$phospho_r, c("Gene HGNC symbol","Gene description", "Site position", "2h", "6h", "12h", "18h", "24h"))
  return(table_data)
}

