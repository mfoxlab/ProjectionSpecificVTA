workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/NAC/iRegulon"
NACworkdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/NAC"
setwd(workdir)

library(dplyr)
library(readxl)
library(stringr) # For extracting numbers from sheet names


NAC_gene_summary <- read.table(file.path(NACworkdir, "NAC_gene_summary.txt"), sep = '\t', header = TRUE, quote = '')

excel_file <- "NAC iRegulon Masterlist.xlsx"
sheets <- excel_sheets(excel_file)



# Loop through each sheet
for (sheet in sheets) {
  
  # Read the current sheet
  sheet_data <- read_excel(excel_file, sheet = sheet)
  
  # Extract the sixth column
  sixth_col <- sheet_data[[6]]
  
  # Split the column data by commas, unlist, and get unique values
  split_tfs <- unique(unlist(strsplit(sixth_col, ",")))
  
  # Clean the data by removing empty strings and unwanted entries
  cleaned_tfs <- split_tfs[split_tfs != "" & split_tfs != "Transcription factor"]
  
  # Convert to a data frame with the column name 'mgi_symbol'
  cleaned_tfs_df <- as.data.frame(cleaned_tfs) %>%
    rename(mgi_symbol = cleaned_tfs) %>%
    na.omit()
  
  # Join with the NAC_gene_summary by 'mgi_symbol'
  joinedNAC <- inner_join(cleaned_tfs_df, NAC_gene_summary, by = 'mgi_symbol')
  
  # Extract the number from the sheet name (e.g., 'module5' becomes '5')
  module_number <- str_extract(sheet, "\\d+")
  
  # Construct the file name using the extracted number
  output_file <- paste0("NAC_M", module_number, "_TF_DE.txt")
  
  # Write the resulting data to a file
  write.table(joinedNAC, file = output_file, row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
  
  # Print the file name for reference
  cat("Exported:", output_file, "\n")
}
