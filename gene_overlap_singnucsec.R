# Addiction Biology Paper for cluster info:
# https://onlinelibrary.wiley.com/doi/full/10.1111/adb.13403

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/SING NUC SEQ OVERLAP"

setwd(workdir)


library(dplyr)
library(data.table)


load("~/WGCNA/LibraDE_MouseRatclusters.rda")

# CODE FOR OVERLAP ----------------------------------------------


# inner join the DE gene with my projection specific genes. this will show me specific types of nuclei that these genes are. Via clusters


# METHOD 1. AMY_gene_summary. Has fewer entries so I am just skipping this because idk why it's been filtered out.

# AMYgenes_genesummary <- read.table("AMY/AMY_gene_summary.txt", sep = '\t', header = TRUE, quote = '')



# METHOD 2 - RRHO original data. WORKS.
# USE THIS BECAUSE 'Here is the dataframe that contains all the DEGs regardless of p value for RRHO'

load("/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/RRHO/DEoutputNOP.Rdata")

# Initialize an empty list to store extracted columns
data_list <- list()

# Iterate over each REGION (df_name)
for (df_name in names(DEoutp)) {
  # Get the dataframe
  df <- DEoutp[[df_name]]
  
  # Define dynamic column names
  p_value_col <- paste0("P.Value_", df_name, "_FENT")
  logFC_col <- paste0("logFC_", df_name, "_FENT")
  
  print(p_value_col)
  
  # Select only the P-Value and logFC columns
  region_data <- df[, c(p_value_col, logFC_col), drop = FALSE]
  
  # Rename columns to indicate their region
  colnames(region_data) <- c(paste0("P.Value_", df_name, "_FENT"), paste0("logFC_", df_name, "_FENT"))
  
  print(colnames(region_data))
  
  # Store in list for merging later
  data_list[[df_name]] <- region_data
}

# Extract first 4 columns from the first dataframe in DEoutp (assumes all have the same structure)
base_df <- DEoutp[[names(DEoutp)[1]]][, c("ensembl_gene_id", "mgi_symbol", "gene_biotype", "description")]

# Combine all extracted data into one final dataframe
all_region_genes <- cbind(base_df, do.call(cbind, data_list))

# Export the final dataframe as a tab-delimited text file
write.table(all_region_genes, "all_region_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#create a column called gene for easy inner_join of the data sets because in the original subset its called mgi_symbol 
DE <- DE %>%
  rename_at("gene", ~"mgi_symbol")

#join the two datasets by their common values in gene 
overlapped_genes <-inner_join(all_region_genes, DE, by = "mgi_symbol")


clusternum <- list(16, 18, 13, 9, 5, 11, 15, 20, 24, 19)

# 16 - DA
# 18 - GLUGABADA
# 13 - GLUTAMATE
#  9 - GABA
#  5 - GLU-GABA ISH 
# 11 GLUTAMATE & VMAT2
# 15 - MOSTLY GABA
# 20 - GLUTAMATE
# 24 - GABA
# 19 - GABADA

for (num in clusternum) {
  
  # name for each filtered data frame
  currentcluster <- paste0("totalgenescluster", num) 
  
  # make the filtered data frame of the cluster number
  filtered_df <- data.frame(overlapped_genes %>%
    filter(cell_type == num,
           p_val < 0.05))
  
  # export to .txt files
  write.table(assign(currentcluster, filtered_df), paste0(currentcluster, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # PERMUTATION 2: FILTER BY PROJECTION
  
  # AMY
  currentclusterAMY <- paste0(currentcluster, "_AMY")
  filtered_df_proj <- filtered_df[,c(-5,-6,-9,-10,-11,-12)]
  print(filtered_df_proj)
  #write.table(assign(currentclusterAMY, filtered_df_proj), paste0(currentclusterAMY, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterAMY, filtered_df_proj))
  
  
  # PFC
  currentclusterPFC <- paste0(currentcluster, "_PFC")
  filtered_df_proj <- filtered_df[,c(-7:-12)]
  print(filtered_df_proj)
  #write.table(assign(currentclusterPFC, filtered_df_proj), paste0(currentclusterPFC, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterPFC, filtered_df_proj))
  
  
  # NAC
  currentclusterNAC <- paste0(currentcluster, "_NAC")
  filtered_df_proj <- filtered_df[,c(-5:-10)]
  print(filtered_df_proj)
  #write.table(assign(currentclusterNAC, filtered_df_proj), paste0(currentclusterNAC, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterNAC, filtered_df_proj))
  
  # INPUT
  currentclusterINPUT <- paste0(currentcluster, "_INPUT")
  filtered_df_proj <- filtered_df[,c(-5:-10)]
  print(filtered_df_proj)
  #write.table(assign(currentclusterINPUT, filtered_df_proj), paste0(currentclusterINPUT, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterINPUT, filtered_df_proj))
  
  
}

overlapped_genes %>%
  filter(cell_type == 19) # no results.


'
# NOW OVERLAP THE TWO BY GENE NAME
AMYgenes <- read.table("AMY_filtered.txt", sep = '\t', header = TRUE, quote = '')

#create a column called gene for easy inner_join of the data sets because in the original subset its called mgi_symbol 
AMYgenes$gene <- AMYgenes$mgi_symbol

#join the two datasets by their common values in gene 
AMYclustergenes <-inner_join(AMYgenes, DE, by = "gene")

'


save.image(file = paste0(workdir, "/gene_overlap_singnucseq_workspace.Rdata"))











