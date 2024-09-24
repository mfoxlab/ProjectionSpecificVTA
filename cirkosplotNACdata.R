# CREATING A FILE FOR CIRKOS PLOTS

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/NAC"
setwd(workdir)

# Load the workspace
load(file = paste0(workdir,"/cirkosplotNACdataworkspace.Rdata"))

library(dplyr)

copieddata <- read.delim("clipboard", header = TRUE)


NAC_gene_summary <- read.table("NAC_gene_summary.txt", sep = '\t', header = TRUE, quote = '')

joineddf <- inner_join(copieddata, NAC_gene_summary, by = 'ensembl_gene_id') 

joineddf <- joineddf %>%
  select('ModuleEigengene', 'mgi_symbol', 'ensembl_gene_id', 'logFC_NAC_FENT', 'logFC_NAC_F_FENT', 'logFC_NAC_M_FENT' , 'P.Value_NAC_FENT', 'P.Value_NAC_F_FENT', 'P.Value_NAC_M_FENT') 


# EXPORT AS AN EXCEL FILE
write.table(joineddf, file = "NAC_cirkosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)



# Save the workspace
save.image(file = paste0(workdir, "/cirkosplotNACdataworkspace.Rdata"))
