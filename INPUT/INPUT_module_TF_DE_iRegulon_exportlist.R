workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/INPUT/iRegulon"
INPUTworkdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/INPUT"
setwd(workdir)

library(dplyr)

INPUT_gene_summary <- read.table(file.path(INPUTworkdir, "INPUT_gene_summary.txt"), sep = '\t', header = TRUE, quote = '')
copied_tfs <- read.delim("clipboard", header = FALSE, sep = ",")

copied_tfs <- unique(unlist(copied_tfs))
cleaned_tfs <- unique_tfs[unique_tfs != "" & unique_tfs != "Transcription factor"]
cleaned_tfs <- as.data.frame(cleaned_tfs) %>%
  rename(mgi_symbol = cleaned_tfs)

joinedINPUT <- inner_join(cleaned_tfs, INPUT_gene_summary, by = 'mgi_symbol') 

write.table(joinedINPUT, file = "INPUT_M23_TF_DE.txt", row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
