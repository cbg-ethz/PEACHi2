###################
# Load the filtered data
##################

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

extract_phospho_df = function (ph_table) {
  phospho_names = as.vector(sapply(unlist(ph_table[,1]), function(x) strsplit(as.character(x),"_")[[1]][1]))
  phospho_sites = as.vector(sapply(unlist(ph_table[,1]), function(x) strsplit(as.character(x),"_")[[1]][2]))
  new_table = as.data.frame(cbind(phospho_names, phospho_sites, ph_table[,2:ncol(ph_table)]))
  colnames(new_table) = c("Gene_name", "Site_position", "ratio_2h", "ratio_6h", "ratio_12h", "ratio_18h", "ratio_24h")
  return(new_table)
}

# initialization
rna = prot = phospho = NULL

#RNA
rna$All=as.data.frame(read.table("data_new/initial_ratios/ratios_RNA_final.txt", sep="\t", row.names=NULL, header=TRUE, as.is=TRUE))
prot$All = as.data.frame(read.table("data_new/initial_ratios/ratios_proteins_filtered.txt", sep="\t", header=TRUE))
phospho$All = as.data.frame(read.table("data_new/initial_ratios/ratios_phospho_normalized.txt", sep="\t", row.names=NULL, header=TRUE, as.is=TRUE, stringsAsFactors = FALSE))

# cutoffs
rna$cutoff = read.table("data_new/DE/DE_RNA_cutoff_2.txt")
prot$cutoff = read.table("data_new/DE/DE_prot_cutoff_2.txt")
phospho$cutoff = read.table("data_new/DE/DE_phospho_cutoff_2.txt")

# counts
rna$counts = as.data.frame(read.table("data_new/initial_ratios/counts_RNA_normalized.txt", sep="\t", row.names=NULL, header=TRUE, as.is=TRUE))
prot$counts = as.data.frame(read.table("data_new/initial_ratios/ratios_counts_proteins_filtered.txt", sep=" ", row.names=NULL, header=TRUE, as.is=TRUE))
prot$counts = prot$counts[,-c(2:6)]

# annotations
gene_annotations = read.table("data_new/gene_annotation_table.txt", header=TRUE, row.names=NULL, sep="\t", quote = "", as.is=TRUE, stringsAsFactors = FALSE)

# interesting groups
rna$Up = read.table("data_new/clusters_ratios/rna_Up.txt", row.names=NULL)
rna$Down = read.table("data_new/clusters_ratios/rna_Down.txt", row.names=NULL)
rna$Not = rna$All[which(rna$All[,1] %in% setdiff(rna$All[,1],union(rna$Up[,1],rna$Down[,1]))),]

prot$Up = read.table("data_new/clusters_ratios/prot_Up.txt", row.names=NULL)
prot$Down = read.table("data_new/clusters_ratios/prot_Down.txt", row.names=NULL)
prot$Not = prot$All[which(prot$All[,1] %in% setdiff(prot$All[,1],union(prot$Up[,1],prot$Down[,1]))),]

phospho$Up = extract_phospho_df(read.table("data_new/clusters_ratios/phospho_Up.txt", row.names=NULL))
phospho$Down = extract_phospho_df(read.table("data_new/clusters_ratios/phospho_Down.txt", row.names=NULL))
p_up_down = rbind(phospho$Up, phospho$Down)
phospho$Not = phospho$All[which(sapply(paste(phospho$All[,1], phospho$All[,2], sep="_"), 
                                          function(x) !(x %in% paste(p_up_down[,1], p_up_down[,2], sep="_")))),]