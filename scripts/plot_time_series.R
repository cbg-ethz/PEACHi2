# plots the time series from a cluster as well as their mean and standard deviation
library(zoo)
plot_time_series = function(time_series_matrix, cutoff_matrix, time_points, yl, y_label, x_label, titl, curve_color, s_time_series) {
  curve_color_a = adjustcolor(curve_color, alpha.f = 0.8)
  if (nrow(time_series_matrix) > 0) {
    # plot the time series
    for( i in 1:nrow(time_series_matrix) ) {
      plot(time_points, time_series_matrix[i,], ylim=yl, xlim = c(0,30), axes = FALSE, xlab="", ylab="", main="", col=curve_color_a, type="l", pch=20,lwd=3)
      par(new=TRUE)
      n_na = which(!is.na(time_series_matrix[i,]))
      if(sum(n_na < 3)>0 & sum(n_na > 3)>0 & length(n_na) < 5 & (1 %in% n_na)) {
        plot(time_points, na.spline(c(time_series_matrix[i,])), ylim=yl, xlim = c(0,30), axes = FALSE, xlab="", ylab="", main="", col=curve_color_a, type="l", lty=3, lwd=3)
        par(new=TRUE)
      }
    }
    
    # plot the DE cutoffs
    plot(time_points, cutoff_matrix[1,], ylim=yl, xlim = c(0,30), axes=FALSE, xlab = "", ylab = "", main = "", col = adjustcolor("darkgrey",alpha.f=0.7), type="l", lty = 2, lwd = 4)
    par(new=TRUE)
    plot(time_points, cutoff_matrix[2,], ylim=yl, xlim = c(0,30), axes=FALSE, xlab = "", ylab = "", main = "", col = adjustcolor("darkgrey",alpha.f=0.7), type="l", lty = 2, lwd = 4)
    par(new=TRUE)
    
    # plot the selected genes
    for( i in 1:nrow(s_time_series) ) {
      plot(time_points, s_time_series[i,], ylim=yl, xlim = c(0,30), axes = FALSE, xlab="", ylab="", main="", col="black", type="l", pch=20,lwd=2)
      par(new=TRUE)
      n_na = which(!is.na(s_time_series[i,]))
      if(sum(n_na < 3)>0 & sum(n_na > 3)>0 & length(n_na) < 5 & (1 %in% n_na)) {
        plot(time_points, na.spline(c(s_time_series[i,])), ylim=yl, xlim = c(0,30), axes = FALSE, xlab="", ylab="", main="", col="black", type="l", lty=3, lwd=2)
        par(new=TRUE)
      }
    }
    
    # finalize plot details (axes labels, title)
    axis(side = 1, at = c(0,2,6,12,18,24,30),cex.axis=2.5,cex.lab=2.5)
    axis(side = 2,cex.axis=2.5,cex.lab=2.5)
    title(main=titl, xlab = x_label, ylab = y_label,cex.lab=2.5,cex.main=2.5)
    par(new=FALSE)
  }
}

print_error_message = function(message) {
  plot(0,0,xlim = c(0,30), ylim=c(0,6), axes = FALSE, ylab = "", xlab = "", type = "l")
  text(x = 10, y=2, labels = message, cex = 2)
}

plot_trio_given_list = function(rna, prot, phospho, rna_c, prot_c, phospho_c, gene_list, color_array, selected_rows) {
  # display an error message if the requested behavior combination has no genes
  if (nrow(gene_list) == 0) {
    print_error_message("There are no genes with the selected patterns.")
    return()
  }
  
  # extract the genes time series
  if (sum(gene_list == "all")==1) {
    cluster_rna_ratios = rna[,2:ncol(rna)]
    cluster_prot_ratios = prot[,2:ncol(prot)]
    cluster_phospho_ratios = phospho[,3:ncol(phospho)]
  } else {
    cluster_rna_ratios = rna[which(toupper(rna[,1]) %in% toupper(gene_list[,1])),2:ncol(rna)]
    cluster_prot_ratios = prot[which(toupper(prot[,1]) %in% toupper(gene_list[,1])),2:ncol(prot)]
    cluster_phospho_ratios = phospho[which(toupper(phospho[,1]) %in% toupper(gene_list[,1])),3:ncol(phospho)]
  }
  
  # extract the selected genes time series
  selected_genes = c(selected_rows$rna_selected[,1], selected_rows$prot_selected[,1], selected_rows$ph_selected[,1])
  selected_rna_ratios = rna[which(rna[,1] %in% selected_genes),2:ncol(rna)]
  selected_prot_ratios = prot[which(prot[,1] %in% selected_genes),2:ncol(prot)]
  selected_phospho_ratios = phospho[which(phospho[,1] %in% selected_genes),3:ncol(phospho)]
  
  # display error message if the none of the requested genes are found
  # if (nrow(cluster_rna_ratios) == 0 & nrow(cluster_prot_ratios) == 0 & nrow(cluster_phospho_ratios) == 0) {
  #   print_error_message("There were no genes with the selected behavior pattern")
  #   return()
  # }
  
  par(mfrow=c(1,3), mar=c(5,6,4,2)+0.1) 
  range_plot = c(-4,6) #range(cluster_rna_ratios,cluster_prot_ratios,cluster_phospho_ratios,rna_c, prot_c, phospho_c, na.rm = TRUE)
  if(nrow(cluster_rna_ratios)>0){
    plot_time_series(cluster_rna_ratios, rna_c, c(6,12,18,24,30), range_plot, "log2(HIV/Mock)", "", "RNA", color_array[1], selected_rna_ratios)
  }
  if(nrow(cluster_prot_ratios)>0){
    plot_time_series(cluster_prot_ratios, prot_c, c(2,6,12,18,24), range_plot, "", "time(h)", "Protein", color_array[2], selected_prot_ratios)
  }
  if(nrow(cluster_phospho_ratios)>0){
    plot_time_series(cluster_phospho_ratios, phospho_c, c(2,6,12,18,24), range_plot, "", "", "Phosphoprotein", color_array[3], selected_phospho_ratios)
  }
}




