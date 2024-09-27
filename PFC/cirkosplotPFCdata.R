# CREATING A FILE FOR CIRKOS PLOTS

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/PFC"
setwd(workdir)

# Load the workspace
load(file = paste0(workdir,"/cirkosplotPFCdataworkspace.Rdata"))

library(dplyr)

copieddata <- read.delim("clipboard", header = TRUE)


PFC_gene_summary <- read.table("PFC_gene_summary.txt", sep = '\t', header = TRUE, quote = '')

joineddf <- inner_join(copieddata, PFC_gene_summary, by = 'ensembl_gene_id') 

joineddf <- joineddf %>%
  select('ModuleEigengene', 'mgi_symbol', 'ensembl_gene_id', 'logFC_PFC_FENT', 'logFC_PFC_F_FENT', 'logFC_PFC_M_FENT' , 'P.Value_PFC_FENT', 'P.Value_PFC_F_FENT', 'P.Value_PFC_M_FENT') 


# EXPORT AS AN EXCEL FILE
write.table(joineddf, file = "PFC_cirkosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)



# Save the workspace
save.image(file = paste0(workdir, "/cirkosplotPFCdataworkspace.Rdata"))
