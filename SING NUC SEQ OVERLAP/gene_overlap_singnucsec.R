# Addiction Biology Paper for cluster info:
# https://onlinelibrary.wiley.com/doi/full/10.1111/adb.13403

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/SING NUC SEQ OVERLAP"

setwd(workdir)

# install libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)


load(file = paste0(workdir, "/gene_overlap_singnucseq_workspace.Rdata"))
load("~/WGCNA/LibraDE_MouseRatclusters.rda")

length(unique(DE$mgi_symbol)) #49147 unique genes.


# CODE FOR OVERLAP ----------------------------------------------
# inner join the DE gene with my projection specific genes. this will show me specific types of nuclei that these genes are. Via clusters


# METHOD 1. AMY_gene_summary. Has fewer entries so I am just skipping this because idk why it's been filtered out.

# AMYgenes_genesummary <- read.table("AMY/AMY_gene_summary.txt", sep = '\t', header = TRUE, quote = '')



# METHOD 2 - RRHO original data. WORKS.
# USE THIS BECAUSE 'Here is the dataframe that contains all the DEGs regardless of p value for RRHO'

'''
test <- (DEoutp$PFC$P.Value_AMY_PFC_SAL) < 0.05
summary(test)
'''

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
overlapped_genes <- inner_join(all_region_genes, DE, by = "mgi_symbol") %>%
  filter(p_val < 0.05)

clusternum <- list(16, 18, 13, 9, 5, 11, 15, 20, 24, 19)

# 16 - DA
# 18 - GLUGABADA
# 13 - GLUTAMATE
#  9 - GABA
#  5 - GLU-GABA ISH 
# 11 - GLUTAMATE & VMAT2
# 15 - MOSTLY GABA
# 20 - GLUTAMATE
# 24 - GABA
# 19 - GABADA

for (num in clusternum) {
  
  # name for each filtered data frame
  currentcluster <- paste0("totalgenescluster", num) 
  
  # make the filtered data frame of the cluster number
  filtered_df <- data.frame(overlapped_genes %>%
    filter(cell_type == num))
        
  
  # export to .txt files
  write.table(assign(currentcluster, filtered_df), paste0(currentcluster, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # PERMUTATION 2: FILTER BY PROJECTION
  
  # AMY
  currentclusterAMY <- paste0(currentcluster, "_AMY")
  filtered_df_proj <- filtered_df[,c(-5,-6,-9,-10,-11,-12)]
  #write.table(assign(currentclusterAMY, filtered_df_proj), paste0(currentclusterAMY, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterAMY, filtered_df_proj))
  
  
  # PFC
  currentclusterPFC <- paste0(currentcluster, "_PFC")
  filtered_df_proj <- filtered_df[,c(-7:-12)]
  #write.table(assign(currentclusterPFC, filtered_df_proj), paste0(currentclusterPFC, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterPFC, filtered_df_proj))
  
  
  # NAC
  currentclusterNAC <- paste0(currentcluster, "_NAC")
  filtered_df_proj <- filtered_df[,c(-5:-10)]
  #write.table(assign(currentclusterNAC, filtered_df_proj), paste0(currentclusterNAC, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterNAC, filtered_df_proj))
  
  # INPUT
  currentclusterINPUT <- paste0(currentcluster, "_INPUT")
  filtered_df_proj <- filtered_df[,c(-5:-10)]
  #write.table(assign(currentclusterINPUT, filtered_df_proj), paste0(currentclusterINPUT, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.frame(assign(currentclusterINPUT, filtered_df_proj))
  
  
}

overlapped_genes %>%
  filter(cell_type == 19) # no results.


overlapped_genes <- overlapped_genes %>%
  mutate(
    cell_type, cell_type = recode(cell_type,
    '16' = '16: DA',
    '18' =  '18: GLUGABADA',
    '13' = '13: GLU',
    '9' =  '9: GABA',
    '5' = '5: GLUGABA ISH',
    '11' = '11: GLU + VMAT2',
    '15' =  '15: MOSTLY GABA',
    '20' =  '20: GLU',
    '24' = '24: GABA',
    '19'  = '19: GABADA',
    '14' = '14: GABA'
  ))


'
# NOW OVERLAP THE TWO BY GENE NAME
AMYgenes <- read.table("AMY_filtered.txt", sep = "\t", header = TRUE, quote = '')

#create a column called gene for easy inner_join of the data sets because in the original subset its called mgi_symbol 
AMYgenes$gene <- AMYgenes$mgi_symbol

#join the two datasets by their common values in gene 
AMYclustergenes <-inner_join(AMYgenes, DE, by = "gene")

'


# INITIAL DATA VISUALIZATION ----------------------------------------------
# 1: bar graph - total cluster makeup.

graph1_total_cluster_counts <- ggplot(data = overlapped_genes, aes(x = reorder(cell_type, -table(cell_type)[cell_type]), fill = cell_type)) +
  geom_bar() +
  labs(title = "Total Cluster Makeup (p < 0.05)", x = "Cell Type", y = "Nuclei Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") + # Adjust colors for better visualization
  geom_text(stat = 'count', aes(label = after_stat(count)), size = 4, position = position_stack(vjust = 0.5))

graph1_total_cluster_counts


overlapped_genes %>%
  count(cell_type) %>% # gaba is the majority. Matches the subtables.
  arrange(desc(n))

all_region_genes %>%
  count(gene_biotype) # the diversity of the types of genes regardless of significance


length(unique(overlapped_genes$mgi_symbol)) #5584 unique genes. Does this mean that I am plotting multiple genes? (No, I am plotting nuclei!)

# make a new table and filter by region and count the number of sig genes per region???


# 16 - DA
# 18 - GLUGABADA
# 13 - GLUTAMATE
#  9 - GABA
#  5 - GLU-GABA ISH 
# 11 - GLUTAMATE & VMAT2
# 15 - MOSTLY GABA
# 20 - GLUTAMATE
# 24 - GABA
# 19 - GABADA

# 2: stacked bar graph - cluster makeup by region

cluster_makeup_PFC <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  count(cell_type) %>%
  remove_rownames %>%
  column_to_rownames(var = "cell_type") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region")
cluster_makeup_PFC$region <- "PFC"


cluster_makeup_NAC <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  count(cell_type) %>%
  remove_rownames %>%
  column_to_rownames(var = "cell_type") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region") 
cluster_makeup_NAC$region <- "NAC"


cluster_makeup_AMY <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  count(cell_type) %>%
  remove_rownames %>%
  column_to_rownames(var = "cell_type") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region") 
cluster_makeup_AMY$region <- "AMY"


cluster_makeup_INPUT <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  count(cell_type) %>%
  remove_rownames %>%
  column_to_rownames(var = "cell_type") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region") 
cluster_makeup_INPUT$region <- "INPUT"


# combine
cluster_makeup_REGION <- bind_rows(cluster_makeup_PFC, cluster_makeup_NAC, cluster_makeup_AMY, cluster_makeup_INPUT) %>%
  pivot_longer(cols = -region, names_to = "cell cluster", values_to = "cell count")

cluster_makeup_counts_REGION_order <- cluster_makeup_REGION %>%
  group_by(region) %>%
  summarize('total region count' = sum(`cell count`)) %>%
  arrange(desc(`total region count`)) %>%
  pull(region)

# graph
graph2_cluster_makeup_REGION <- ggplot(cluster_makeup_REGION, aes(x = factor(region, levels = cluster_makeup_counts_REGION_order), y = `cell count`, fill = `cell cluster`)) +
  geom_bar(stat = "identity") +
  labs(title = "Cluster Makeup by Region", x = "Region", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Adjust colors for better visualization 
  geom_text(aes(label = `cell count`), size = 4, position = position_stack(vjust = 0.5))

graph2_cluster_makeup_REGION


# 2.1: now do the same but with decimal ratios 

cluster_makeup_REGION_decimalvalues <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05 |
        NAC.P.Value_NAC_FENT < 0.05 |
          AMY.P.Value_AMY_FENT < 0.05 |
          INPUT.P.Value_INPUT_FENT < 0.05) %>%
  group_by(cell_type) %>%
  count(cell_type) %>%
  rename_at('cell_type', ~'cell cluster')

cluster_makeup_REGION <- cluster_makeup_REGION %>%
  inner_join(cluster_makeup_REGION_decimalvalues, by = 'cell cluster') %>%
  mutate(
    `cell count` = as.numeric(`cell count`),
    n = as.numeric(n),
    ratio = `cell count` / n) %>%
  mutate(ratio = round(ratio, 2))


# graph
graph2.1_cluster_makeup_REGION_decimalvalues <- ggplot(cluster_makeup_REGION, aes(x = factor(region, levels = cluster_makeup_counts_REGION_order), y = ratio, fill = `cell cluster`)) +
  geom_bar(stat = "identity") +
  labs(title = "Cluster Makeup by Region", x = "Region", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Adjust colors for better visualization 
  geom_text(aes(label = ratio), size = 4, position = position_stack(vjust = 0.5))

graph2.1_cluster_makeup_REGION_decimalvalues
  

# 3: stacked bar graph - region makeup with subclusters

REGION_counts <- overlapped_genes %>%
  group_by(cell_type) %>%
  summarise(
    AMY = sum(AMY.P.Value_AMY_FENT < 0.05, na.rm = TRUE),
    PFC = sum(PFC.P.Value_PFC_FENT < 0.05, na.rm = TRUE),
    NAC = sum(NAC.P.Value_NAC_FENT < 0.05, na.rm = TRUE),
    INPUT = sum(INPUT.P.Value_INPUT_FENT < 0.05, na.rm = TRUE)
  ) %>%
  rename(cluster = cell_type) %>%
  as.data.frame() %>%
  pivot_longer(cols = c(AMY, PFC, NAC, INPUT),
               names_to = 'region',
               values_to = 'cell count')

REGION_counts_order <- REGION_counts %>%
  group_by(cluster) %>%
  summarise(total_count = sum(`cell count`)) %>%
  arrange(desc(total_count)) %>%
  pull(cluster)

# graph
graph3_REGION_cluster_makeup <- ggplot(REGION_counts, aes(x = factor(cluster, levels = REGION_counts_order), y = `cell count`, fill = `region`)) +
  geom_bar(stat = "identity") +
  labs(title = "Region Makeup with Subclusters", x = "Region", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Adjust colors for better visualization 
  geom_text(aes(label = `cell count`), size = 4, position = position_stack(vjust = 0.5))

graph3_REGION_cluster_makeup

### MAKING UP-DOWN COLORS

UPDOWN_fillcolors <- c("UP" = "deeppink2", "DOWN" = "deepskyblue3")


# 4: dodged bar graph UP-DOWN VERSION of VTA -> REGION cluster makeup

cluster_makeup_UPDOWN <- overlapped_genes %>%
  group_by(cell_type) 

cluster_makeup_UPDOWN$avg_logFC <- ifelse(cluster_makeup_UPDOWN$avg_logFC < 0, "DOWN", "UP")

cluster_makeup_UPDOWN <- cluster_makeup_UPDOWN %>%
  count(avg_logFC)

# graph
graph4_cluster_makeup_UPDOWN <- ggplot(cluster_makeup_UPDOWN, aes(x = cell_type, y = n, fill = avg_logFC)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "Log Fold Change (Up/Downregulated)", x = "Cell Type", y = "Total Nuclei Count (NO p < 0.05)") +
  theme_minimal() +
  scale_fill_manual(values = UPDOWN_fillcolors) +  
  geom_text(aes(label = n), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph4_cluster_makeup_UPDOWN


# 5: dodged bar graph UP-DOWN VERSION of VTA -> REGION makeup

REGION_counts_UPDOWN <- overlapped_genes

REGION_counts_UPDOWN$avg_logFC <- ifelse(REGION_counts_UPDOWN$avg_logFC < 0, "DOWN", "UP")

REGION_counts_UPDOWN <- REGION_counts_UPDOWN %>%
  group_by(avg_logFC) %>%
  summarise(
    AMY = sum(AMY.P.Value_AMY_FENT < 0.05, na.rm = TRUE),
    PFC = sum(PFC.P.Value_PFC_FENT < 0.05, na.rm = TRUE),
    NAC = sum(NAC.P.Value_NAC_FENT < 0.05, na.rm = TRUE),
    INPUT = sum(INPUT.P.Value_INPUT_FENT < 0.05, na.rm = TRUE)
  ) %>%
  as.data.frame() %>%
  pivot_longer(cols = c(AMY, PFC, NAC, INPUT),
               names_to = 'region',
               values_to = 'cell count')

# graph
graph5_REGION_counts_UPDOWN <- ggplot(REGION_counts_UPDOWN, aes(x = factor(region), y = `cell count`, fill = `avg_logFC`)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "Log Fold Change (Up/Downregulated) by Region", x = "Region", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_manual(values = UPDOWN_fillcolors) +  
  geom_text(aes(label = `cell count`), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph5_REGION_counts_UPDOWN


# 6: region-specific cluster UP-DOWN counts

alluvial_df <- overlapped_genes %>%
  mutate(
    PFC_reg = ifelse(PFC.logFC_PFC_FENT < 0, "DOWN", "UP"),
    AMY_reg = ifelse(AMY.logFC_AMY_FENT < 0, "DOWN", "UP"),
    NAC_reg = ifelse(NAC.logFC_NAC_FENT < 0, "DOWN", "UP"),
    INPUT_reg = ifelse(INPUT.logFC_INPUT_FENT < 0, "DOWN", "UP"),
    avg_logFC_reg = ifelse(avg_logFC < 0, "DOWN", "UP")
  )

# 6.1: PFC-specific cluster UP-DOWN counts
alluvial_PFC <- alluvial_df %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  group_by(cell_type) %>%
  count(avg_logFC_reg)

# graph
graph6.1_PFC_cluster_counts_UPDOWN <- ggplot(alluvial_PFC, aes(x = cell_type, y = n, fill = avg_logFC_reg)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "PFC Log Fold Change (Up/Downregulated)", x = "Cell Type", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_manual(values = UPDOWN_fillcolors) +  
  geom_text(aes(label = n), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph6.1_PFC_cluster_counts_UPDOWN

# 6.2: AMY-specific cluster UP-DOWN counts
alluvial_AMY <- alluvial_df %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  group_by(cell_type) %>%
  count(avg_logFC_reg)

# graph
graph6.2_AMY_cluster_counts_UPDOWN <- ggplot(alluvial_AMY, aes(x = cell_type, y = n, fill = avg_logFC_reg)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "AMY Log Fold Change (Up/Downregulated)", x = "Cell Type", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_manual(values = UPDOWN_fillcolors) +  
  geom_text(aes(label = n), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph6.2_AMY_cluster_counts_UPDOWN

# 6.3: NAC-specific cluster UP-DOWN counts
alluvial_NAC <- alluvial_df %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  group_by(cell_type) %>%
  count(avg_logFC_reg)

# graph
graph6.3_NAC_cluster_counts_UPDOWN <- ggplot(alluvial_NAC, aes(x = cell_type, y = n, fill = avg_logFC_reg)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "NAC Log Fold Change (Up/Downregulated)", x = "Cell Type", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_manual(values = UPDOWN_fillcolors) +  
  geom_text(aes(label = n), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph6.3_NAC_cluster_counts_UPDOWN

# ALLUVIAL PLOTS ----------------------------------------------

# 6.4: INPUT-specific cluster UP-DOWN counts
alluvial_INPUT <- alluvial_df %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  group_by(cell_type) %>%
  count(avg_logFC_reg)

# graph
graph6.4_INPUT_cluster_counts_UPDOWN <- ggplot(alluvial_INPUT, aes(x = cell_type, y = n, fill = avg_logFC_reg)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "INPUT Log Fold Change (Up/Downregulated)", x = "Cell Type", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_manual(values = UPDOWN_fillcolors) +  
  geom_text(aes(label = n), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph6.4_INPUT_cluster_counts_UPDOWN


# 6.5: alluvial plot but I am not filtering the IP inj genes. 

# No. First, make bar plots of the combinations

alluvial_PFC_combos <- alluvial_df %>%
  mutate(IP_IV_PFC_reg_comb = paste(PFC_reg, avg_logFC_reg)) %>%
  count(IP_IV_PFC_reg_comb) %>%
  column_to_rownames(var = "IP_IV_PFC_reg_comb") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region")
  
alluvial_PFC_combos['region'][alluvial_PFC_combos['region'] == 'n'] <- 'PFC'


alluvial_AMY_combos <- alluvial_df %>%
  mutate(IP_IV_AMY_reg_comb = paste(AMY_reg, avg_logFC_reg)) %>%
  count(IP_IV_AMY_reg_comb) %>%
  column_to_rownames(var = "IP_IV_AMY_reg_comb") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region")

alluvial_AMY_combos['region'][alluvial_AMY_combos['region'] == 'n'] <- 'AMY'


alluvial_NAC_combos <- alluvial_df %>%
  mutate(IP_IV_NAC_reg_comb = paste(NAC_reg, avg_logFC_reg)) %>%
  count(IP_IV_NAC_reg_comb) %>%
  column_to_rownames(var = "IP_IV_NAC_reg_comb") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region")

alluvial_NAC_combos['region'][alluvial_NAC_combos['region'] == 'n'] <- 'NAC'


alluvial_INPUT_combos <- alluvial_df %>%
  mutate(IP_IV_INPUT_reg_comb = paste(INPUT_reg, avg_logFC_reg)) %>%
  count(IP_IV_INPUT_reg_comb) %>%
  column_to_rownames(var = "IP_IV_INPUT_reg_comb") %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var = "region")

alluvial_INPUT_combos['region'][alluvial_INPUT_combos['region'] == 'n'] <- 'INPUT'


# combine
alluvial_combined_combos <- bind_rows(alluvial_PFC_combos, alluvial_AMY_combos, alluvial_NAC_combos, alluvial_INPUT_combos) %>%
  pivot_longer(cols = -region, names_to = "UP/DOWN regulated: RNASeq & scSeq", values_to = "count")


library(ggalluvial)

ggplot(data = alluvial_combined_combos,
       aes(axis1 = region, axis2 = `UP/DOWN regulated: RNASeq & scSeq`, y = count)) +
  geom_alluvium(aes(fill = `UP/DOWN regulated: RNASeq & scSeq`)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survery", "Regulated"),
                   expand = c(0.15, 0.05)) +
  theme_void()






# redoing alluvial plos. this time I am separating the counts to be 'Fent IP' and 'Fent IVSA'

alluvial_PFC_combos_v2 <- alluvial_df %>%
  count(PFC_reg, avg_logFC_reg) %>%
  mutate(color_comb = paste(PFC_reg, avg_logFC_reg)) %>%
  rename("IP FENT" = PFC_reg,
         "IVSA FENT" = avg_logFC_reg,
         "count" = n) %>%
  mutate(region = 'PFC')

graph7.1_alluvial_PFC_combos_v2 <- ggplot(data = alluvial_PFC_combos_v2,
       aes(axis1 = `IP FENT`, axis2 = `IVSA FENT`, y = count)) +
  geom_alluvium(aes(fill = `color_comb`)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survery", "Regulated"),
                   expand = c(0.15, 0.05)) +
  theme_void()

graph7.1_alluvial_PFC_combos_v2





alluvial_AMY_combos_v2 <- alluvial_df %>%
  count(AMY_reg, avg_logFC_reg) %>%
  mutate(color_comb = paste(AMY_reg, avg_logFC_reg)) %>%
  rename("IP FENT" = AMY_reg,
         "IVSA FENT" = avg_logFC_reg,
         "count" = n) %>%
  mutate(region = 'AMY')

graph7.2_alluvial_AMY_combos_v2 <- ggplot(data = alluvial_AMY_combos_v2,
                                          aes(axis1 = `IP FENT`, axis2 = `IVSA FENT`, y = count)) +
  geom_alluvium(aes(fill = `color_comb`)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survery", "Regulated"),
                   expand = c(0.15, 0.05)) +
  theme_void()

graph7.2_alluvial_AMY_combos_v2




alluvial_NAC_combos_v2 <- alluvial_df %>%
  count(NAC_reg, avg_logFC_reg) %>%
  mutate(color_comb = paste(NAC_reg, avg_logFC_reg)) %>%
  rename("IP FENT" = NAC_reg,
         "IVSA FENT" = avg_logFC_reg,
         "count" = n) %>%
  mutate(region = 'NAC')

graph7.3_alluvial_NAC_combos_v2 <- ggplot(data = alluvial_NAC_combos_v2,
                                          aes(axis1 = `IP FENT`, axis2 = `IVSA FENT`, y = count)) +
  geom_alluvium(aes(fill = `color_comb`)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survery", "Regulated"),
                   expand = c(0.15, 0.05)) +
  theme_void()

graph7.3_alluvial_NAC_combos_v2







alluvial_INPUT_combos_v2 <- alluvial_df %>%
  count(INPUT_reg, avg_logFC_reg) %>%
  mutate(color_comb = paste(INPUT_reg, avg_logFC_reg)) %>%
  rename("IP FENT" = INPUT_reg,
         "IVSA FENT" = avg_logFC_reg,
         "count" = n) %>%
  mutate(region = 'INPUT')

graph7.4_alluvial_INPUT_combos_v2 <- ggplot(data = alluvial_INPUT_combos_v2,
                                          aes(axis1 = `IP FENT`, axis2 = `IVSA FENT`, y = count)) +
  geom_alluvium(aes(fill = `color_comb`)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survery", "Regulated"),
                   expand = c(0.15, 0.05)) +
  theme_void()

graph7.4_alluvial_INPUT_combos_v2


# ok not really the most informative graph. going to make a bar graph:

# graph
graph7_REGION_counts_alluvial_draft <- ggplot(alluvial_combined_combos, aes(x = region, y = `count`, fill = `UP/DOWN regulated: RNASeq & scSeq`)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = 'single' )) +
  labs(title = "Up/Downregulated between IP and IV Fent Genes", x = "Region", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  geom_text(aes(label = `region`), size = 4, position = position_dodge(width = 1), vjust = -0.5)

graph7_REGION_counts_alluvial_draft



graph2_cluster_makeup_REGION <- ggplot(cluster_makeup_REGION, aes(x = factor(region, levels = cluster_makeup_counts_REGION_order), y = `cell count`, fill = `cell cluster`)) +
  geom_bar(stat = "identity") +
  labs(title = "Cluster Makeup by Region", x = "Region", y = "Total Nuclei Count (p < 0.05)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Adjust colors for better visualization 
  geom_text(aes(label = `cell count`), size = 4, position = position_stack(vjust = 0.5))

graph2_cluster_makeup_REGION


# GENE LISTS FOR METASCAPE ANALYSIS ----------------------------------------------

# go terms for different vta projections. 

metapfc <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  distinct(mgi_symbol)
write.table(metapfc, file = "metapfc.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: neurotransmitter transport, collagen chain trimerization, ECM receptor interaction! n = 629

metaamy <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  distinct(mgi_symbol)
write.table(metaamy, file = "metaamy.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: cell junction organization, cell morphogenesis involved in neuron differentiation, renal system development! n = 583

metanac <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  distinct(mgi_symbol)
write.table(metanac, file = "metanac.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: circadian rhythm of gene expression, actin filament organization, regulation of vacuole organization! n = 249

metainput <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  distinct(mgi_symbol)
write.table(metainput, file = "metainput.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: learning/memory, regulation of nervous system development, positive regulation of prostaglandin secretion! n = 243


# filtering these gene lists so that the top 3 cell types are shown only. (nac and input have 4 because the bottom two cell-types are tied for third)

metapfc_topcelltypes <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH' |
         cell_type == '13: GLU'|
         cell_type == '9: GABA') %>%
  distinct(mgi_symbol)
write.table(metapfc_topcelltypes, file = "metapfc_topcelltypes.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: synaptic vesicle cycle, tube morphogenesis, amino biosynthetic process! n = 291

metaamy_topcelltypes <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH' |
           cell_type == '13: GLU'|
           cell_type == '11: GLU + VMAT2') %>%
  distinct(mgi_symbol)
write.table(metaamy_topcelltypes, file = "metaamy_topcelltypes.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: cytoskeleton in muscle cells, cell morphogenesis involved in neuron differentiation, regulation of cellular response to transforming growth factor beta stimulus! n = 258

metanac_topcelltypes <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '20: GLU' |
           cell_type == '5: GLUGABA ISH' |
           cell_type == '13: GLU'|
           cell_type == '11: GLU + VMAT2') %>%
  distinct(mgi_symbol)
write.table(metanac_topcelltypes, file = "metanac_topcelltypes.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: regulation of vacuole organizion, circadian regulation of gene expression, regulation of protein-containing complex disassembly! n = 138

metainput_topcelltypes <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '9: GABA' |
           cell_type == '15: MOSTLY GABA' |
           cell_type == '13: GLU'|
           cell_type == '16: DA') %>%
  distinct(mgi_symbol)
write.table(metainput_topcelltypes, file = "metainput_topcelltypes.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: learning/memory, modulation of chemical synaptic transmission, regulation of nervous system development! n = 130

# filtering even further so that I can see GO terms for each subcluster per projection

# PFC
metapfc_cluster13 <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  distinct(mgi_symbol)
write.table(metapfc_cluster13, file = "metapfc_cluster13.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: morphine addiction, negative regulation of cell-substrate adhesion, PLC-activating GPCR signaling pathway! n = 101

metapfc_cluster9 <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '9: GABA') %>%
  distinct(mgi_symbol)
write.table(metapfc_cluster9, file = "metapfc_cluster9.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: blood vessel development, breast cancer, negative modulation of Notch signaling pathway! n = 96

metapfc_cluster5 <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH') %>%
  distinct(mgi_symbol)
write.table(metapfc_cluster5, file = "metapfc_cluster5.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: synaptic vesicle cycle, neurotransmitter loading into synaptic vesicle, signaling by Receptor Tyrosine Kinases! n = 119

# AMY
metaamy_cluster13 <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  distinct(mgi_symbol)
write.table(metaamy_cluster13, file = "metaamy_cluster13.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: renal system development, regulation of synapse organization, regulation of transforming growth factor beta receptor signaling pathway! n = 88

metaamy_cluster5 <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH') %>%
  distinct(mgi_symbol)
write.table(metaamy_cluster5, file = "metaamy_cluster5.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: neurotransmitter loading into synaptic vesicle, lung development, brain development! n = 100

metaamy_cluster11 <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '11: GLU + VMAT2') %>%
  distinct(mgi_symbol)
write.table(metaamy_cluster11, file = "metaamy_cluster11.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: hormone-mediated signaling pathway, regulation of sodium-ion transmembrane transport, proteoglycan biosynthetic process! n = 86

# NAC
metanac_cluster13 <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  distinct(mgi_symbol)
write.table(metanac_cluster13, file = "metanac_cluster13.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: actin filament organization, puring-containing compound metabolic process! n = 38

metanac_cluster11 <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '11: GLU + VMAT2') %>%
  distinct(mgi_symbol)
write.table(metanac_cluster11, file = "metanac_cluster11.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: negative regulation of glial cell differentiation, metabolism of RNA, cellular response to metal ion! n = 38

metanac_cluster20 <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '20: GLU') %>%
  distinct(mgi_symbol)
write.table(metanac_cluster20, file = "metanac_cluster20.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: DNA repair! n = 38

metanac_cluster5 <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH') %>%
  distinct(mgi_symbol)
write.table(metanac_cluster5, file = "metanac_cluster5.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: regulation of transcription regulatory region DNA binding, cardiac chamber development, regulation of circadian rhythm! n = 37


# INPUT
metainput_cluster13 <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  distinct(mgi_symbol)
write.table(metainput_cluster13, file = "metainput_cluster13.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: cushing syndrome, positive regulation of nervous system development, response to xenobiotic stimulus! n = 41

metainput_cluster9 <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '9: GABA') %>%
  distinct(mgi_symbol)
write.table(metainput_cluster9, file = "metainput_cluster9.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: cytokine-cytokine receptor interaction, regulation of myelination, regulation of nervous system development! n = 40

metainput_cluster15 <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '15: MOSTLY GABA') %>%
  distinct(mgi_symbol)
write.table(metainput_cluster15, file = "metainput_cluster15.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: cell junction organization, regulation of nervous system process, glial cell differentiation! n = 33

metainput_cluster16 <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '16: DA') %>%
  distinct(mgi_symbol)
write.table(metainput_cluster16, file = "metainput_cluster16.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: postsynaptic modulation of chemical synaptic transmission, endochondral ossification, inner ear development! n = 30


# GENE LISTS FOR METASCAPE ANALYSIS - CONCORDANT ----------------------------------------------

# PFC
metapfc_cluster13_concordant <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  filter((PFC.logFC_PFC_FENT < 0 & avg_logFC < 0) |
           PFC.logFC_PFC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metapfc_cluster13_concordant, file = "metapfc_cluster13_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: neural crest differentiation, positive regulation of cell adhesion, blood vessel development! n = 33

metapfc_cluster9_concordant <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '9: GABA') %>%
  filter((PFC.logFC_PFC_FENT < 0 & avg_logFC < 0) |
           PFC.logFC_PFC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metapfc_cluster9_concordant, file = "metapfc_cluster9_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: DNA-regulated transcription initiation! n = 47

metapfc_cluster5_concordant <- overlapped_genes %>%
  filter(PFC.P.Value_PFC_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH') %>%
  filter((PFC.logFC_PFC_FENT < 0 & avg_logFC < 0) |
           PFC.logFC_PFC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metapfc_cluster5_concordant, file = "metapfc_cluster5_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: focal adhesion in P13K Akt mTOR signaling pathway, cardiac ventricle development, purine metabolism! n = 41

# AMY
metaamy_cluster13_concordant <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  filter((AMY.logFC_AMY_FENT < 0 & avg_logFC < 0) |
           AMY.logFC_AMY_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metaamy_cluster13_concordant, file = "metaamy_cluster13_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: regulation of synpase organization, endoplasmic reticulum organization, regulation of blood vessel endothelial cell migration! n = 48

metaamy_cluster5_concordant <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH') %>%
  filter((AMY.logFC_AMY_FENT < 0 & avg_logFC < 0) |
           AMY.logFC_AMY_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metaamy_cluster5_concordant, file = "metaamy_cluster5_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: positive regulation of cell projection organization, establishment of protein localization to membrane, regulation of cellular response to growth factor stimulus! n = 30

metaamy_cluster11_concordant <- overlapped_genes %>%
  filter(AMY.P.Value_AMY_FENT < 0.05) %>%
  filter(cell_type == '11: GLU + VMAT2') %>%
  filter((AMY.logFC_AMY_FENT < 0 & avg_logFC < 0) |
           AMY.logFC_AMY_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metaamy_cluster11_concordant, file = "metaamy_cluster11_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: hormone-mediated signaling pathway, relaxin signaling pathway! n = 35

# NAC
metanac_cluster13_concordant <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  filter((NAC.logFC_NAC_FENT < 0 & avg_logFC < 0) |
           NAC.logFC_NAC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metanac_cluster13_concordant, file = "metanac_cluster13_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: N/A. n = 16

metanac_cluster11_concordant <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '11: GLU + VMAT2') %>%
  filter((NAC.logFC_NAC_FENT < 0 & avg_logFC < 0) |
           NAC.logFC_NAC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metanac_cluster11_concordant, file = "metanac_cluster11_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: protein ubiquitination! n = 18

metanac_cluster20_concordant <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '20: GLU') %>%
  filter((NAC.logFC_NAC_FENT < 0 & avg_logFC < 0) |
           NAC.logFC_NAC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metanac_cluster20_concordant, file = "metanac_cluster20_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: N/A. n = 17

metanac_cluster5_concordant <- overlapped_genes %>%
  filter(NAC.P.Value_NAC_FENT < 0.05) %>%
  filter(cell_type == '5: GLUGABA ISH') %>%
  filter((NAC.logFC_NAC_FENT < 0 & avg_logFC < 0) |
           NAC.logFC_NAC_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metanac_cluster5_concordant, file = "metanac_cluster5_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: N/A. n = 16


# INPUT
metainput_cluster13_concordant <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '13: GLU') %>%
  filter((INPUT.logFC_INPUT_FENT < 0 & avg_logFC < 0) |
           INPUT.logFC_INPUT_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metainput_cluster13_concordant, file = "metainput_cluster13_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: transport of small molecules! n = 13

metainput_cluster9_concordant <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '9: GABA') %>%
  filter((INPUT.logFC_INPUT_FENT < 0 & avg_logFC < 0) |
           INPUT.logFC_INPUT_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metainput_cluster9_concordant, file = "metainput_cluster9_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: cytokine-cytokine receptor interaction, regulation of nervous system development, regulation of hemopoesis! n = 20

metainput_cluster15_concordant <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '15: MOSTLY GABA') %>%
  filter((INPUT.logFC_INPUT_FENT < 0 & avg_logFC < 0) |
           INPUT.logFC_INPUT_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metainput_cluster15_concordant, file = "metainput_cluster15_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: N/A. n = 15

metainput_cluster16_concordant <- overlapped_genes %>%
  filter(INPUT.P.Value_INPUT_FENT < 0.05) %>%
  filter(cell_type == '16: DA') %>%
  filter((INPUT.logFC_INPUT_FENT < 0 & avg_logFC < 0) |
           INPUT.logFC_INPUT_FENT > 0 & avg_logFC > 0) %>%
  distinct(mgi_symbol)
write.table(metainput_cluster16_concordant, file = "metainput_cluster16_concordant.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
# top terms: N/A. n = 12





save.image(file = paste0(workdir, "/gene_overlap_singnucseq_workspace.Rdata"))


