# CREATING A FILE FOR CIRKOS PLOTS

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/INPUT"
setwd(workdir)

# Load the workspace
load(file = paste0(workdir,"/cirkosplotINPUTdataworkspace.Rdata"))

library(dplyr)

copieddata <- read.delim("clipboard", header = TRUE)


INPUT_gene_summary <- read.table("INPUT_gene_summary.txt", sep = '\t', header = TRUE, quote = '')

joineddf <- inner_join(copieddata, INPUT_gene_summary, by = 'ensembl_gene_id') 

joineddf <- joineddf %>%
  select('ModuleEigengene', 'mgi_symbol', 'ensembl_gene_id', 'logFC_INPUT_FENT', 'logFC_INPUT_F_FENT', 'logFC_INPUT_M_FENT' , 'P.Value_INPUT_FENT', 'P.Value_INPUT_F_FENT', 'P.Value_INPUT_M_FENT') 


# EXPORT AS AN EXCEL FILE
write.table(joineddf, file = "INPUT_cirkosdata.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)



# Save the workspace
save.image(file = paste0(workdir, "/cirkosplotINPUTdataworkspace.Rdata"))
