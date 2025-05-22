workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/Annalisa"
setwd(workdir)

load(file = paste0(workdir,"/data_and_figures_workspace.Rdata"))

library(dplyr)
library(ggplot2)

'''
# LOAD IN EACH GENE SUMMARY
PFC_gene_summary <- read.table("/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/PFC/PFC_gene_summary.txt", sep = '\t', header = TRUE, quote = '')
NAC_gene_summary <- read.table("/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/NAC/NAC_gene_summary.txt", sep = '\t', header = TRUE, quote = '')
INPUT_gene_summary <- read.table("/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/INPUT/INPUT_gene_summary.txt", sep = '\t', header = TRUE, quote = '')
'''

PFC_sig_genes <- PFC_gene_summary %>% # filter out all genes that aren't significant per projection
  filter(P.Value_PFC_FENT < 0.05)

NAC_sig_genes <- NAC_gene_summary %>% # filter out all genes that aren't significant per projection
  filter(P.Value_NAC_FENT < 0.05)

INPUT_sig_genes <- INPUT_gene_summary %>% # filter out all genes that aren't significant per projection
  filter(P.Value_INPUT_FENT < 0.05)


# JOINING THE TWO FOR ANNALISA'S VENN DIAGRAM. I initially did by mgi_symbol, but there were gene names that were blank. So now I am doing it by ensembl_gene_id.

# OVERLAP BETWEEN INPUT AND PFC
PFC_INPUT_overlap_genes <- inner_join(PFC_sig_genes, INPUT_sig_genes, by = "ensembl_gene_id") #180 #now 176

# OVERLAP BETWEEN INPUT AND NAC
NAC_INPUT_overlap_genes <- inner_join(NAC_sig_genes, INPUT_sig_genes, by = "ensembl_gene_id") #53 #now 51

# OVERLAP BETWEEN PFC AND NAC
PFC_NAC_overlap_genes <- inner_join(PFC_sig_genes, NAC_sig_genes, by = "ensembl_gene_id") #156. I got a many to many warning message. #now 149



# MAKE 3 TABLES FOR UNION HEATMAPS
# USE 0 TO INDICATE THAT A GENE IS NOT THERE. MAKE IT WHITE OR BLACK

# OUTER/FULL JOIN

# TABLE 1: A TABLE FOR ALL GENES VTA->PFC SIG. INCLUDE THE logFC and P val columns. 
heatmapmerge_PFC_sig_genes <- PFC_sig_genes [,c(1:3, 7, 16)]
heatmapmerge_NAC_sig_genes <- NAC_sig_genes [,c(1:3, 7, 16)]
heatmapmerge_INPUT_sig_genes <- INPUT_sig_genes [,c(1:3, 13, 22)]


heatmapmerged_table <- full_join(heatmapmerge_PFC_sig_genes, heatmapmerge_NAC_sig_genes, by = "ensembl_gene_id") 


# Rename .x columns to _PFC
names(heatmapmerged_table) <- sub("\\.x$", "_PFC", names(heatmapmerged_table))

# Rename .y columns to _NAC
names(heatmapmerged_table) <- sub("\\.y$", "_NAC", names(heatmapmerged_table))

heatmapmerged_table <- full_join(heatmapmerged_table, heatmapmerge_INPUT_sig_genes, by = "ensembl_gene_id") 

# Rename the INPUT columns
heatmapmerged_table <- heatmapmerged_table %>%
  rename_at("mgi_symbol", ~"mgi_symbol_INPUT") %>%
  rename_at("gene_biotype", ~"gene_biotype_INPUT")


heatmapmerged_table <- heatmapmerged_table %>%
  mutate(across(starts_with("logFC"), ~ replace(., is.na(.), 0)))



# MAKING THE TABLE INTO LONG FORMAT

heatmapmerged_table <- heatmapmerged_table %>%
  arrange(logFC_NAC_FENT)

heatmapmerged_table_long <- heatmapmerged_table %>%
  select(ensembl_gene_id, logFC_INPUT_FENT, logFC_PFC_FENT, logFC_NAC_FENT) %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "Region",
    values_to = "logFC"
  ) %>%
  mutate(
    Region = recode(Region,
                    "logFC_INPUT_FENT" = "Total VTA",
                    "logFC_PFC_FENT" = "PFC",
                    "logFC_NAC_FENT" = "NAC"),
    ensembl_gene_id = factor(ensembl_gene_id, levels = heatmapmerged_table$ensembl_gene_id)  # preserve gene order
  )


ggplot(heatmapmerged_table_long, aes(x = ensembl_gene_id, y = Region, fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "deepskyblue3",
    mid = "white",
    high = "deeppink2",
    midpoint = 0,
    limits = c(-2.34, 2.62),
    name = "logFC"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # hides gene labels if too many
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Union Heatmap of Gene Expression"
  )


# table 2: redo but with ALL GENES REGARDLESS OF SIGNIFICANCE
heatmapmerge2_PFC_sig_genes <- PFC_sig_genes [,c(1:3, 7, 16)]
heatmapmerge2_NAC_sig_genes <- NAC_sig_genes [,c(1:3, 7, 16)]
heatmapmerge2_INPUT_sig_genes <- INPUT_sig_genes [,c(1:3, 13, 22)]


heatmapmerged_table <- full_join(heatmapmerge_PFC_sig_genes, heatmapmerge_NAC_sig_genes, by = "ensembl_gene_id") 


# Rename .x columns to _PFC
names(heatmapmerged_table) <- sub("\\.x$", "_PFC", names(heatmapmerged_table))

# Rename .y columns to _NAC
names(heatmapmerged_table) <- sub("\\.y$", "_NAC", names(heatmapmerged_table))

heatmapmerged_table <- full_join(heatmapmerged_table, heatmapmerge_INPUT_sig_genes, by = "ensembl_gene_id") 

# Rename the INPUT columns
heatmapmerged_table <- heatmapmerged_table %>%
  rename_at("mgi_symbol", ~"mgi_symbol_INPUT") %>%
  rename_at("gene_biotype", ~"gene_biotype_INPUT")


heatmapmerged_table <- heatmapmerged_table %>%
  mutate(across(starts_with("logFC"), ~ replace(., is.na(.), 0)))



# MAKING THE TABLE INTO LONG FORMAT

heatmapmerged_table <- heatmapmerged_table %>%
  arrange(logFC_NAC_FENT)

heatmapmerged_table_long <- heatmapmerged_table %>%
  select(ensembl_gene_id, logFC_INPUT_FENT, logFC_PFC_FENT, logFC_NAC_FENT) %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "Region",
    values_to = "logFC"
  ) %>%
  mutate(
    Region = recode(Region,
                    "logFC_INPUT_FENT" = "Total VTA",
                    "logFC_PFC_FENT" = "PFC",
                    "logFC_NAC_FENT" = "NAC"),
    ensembl_gene_id = factor(ensembl_gene_id, levels = heatmapmerged_table$ensembl_gene_id)  # preserve gene order
  )


ggplot(heatmapmerged_table_long, aes(x = ensembl_gene_id, y = Region, fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "deepskyblue3",
    mid = "white",
    high = "deeppink2",
    midpoint = 0,
    limits = c(-2.34, 2.62),
    name = "logFC"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # hides gene labels if too many
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Union Heatmap of Gene Expression"
  )





# Save the workspace
save.image(file = paste0(workdir, "/data_and_figures_workspace.Rdata"))
