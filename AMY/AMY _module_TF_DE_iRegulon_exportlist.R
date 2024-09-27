workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/AMY/iRegulon"
AMYworkdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/AMY"
setwd(workdir)

library(dplyr)
library(readxl)
library(stringr) # For extracting numbers from sheet names


AMY_gene_summary <- read.table(file.path(AMYworkdir, "AMY_gene_summary.txt"), sep = '\t', header = TRUE, quote = '')

excel_file <- "AMY iRegulon Masterlist.xlsx"
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
  
  # Join with the AMY_gene_summary by 'mgi_symbol'
  joinedAMY <- inner_join(cleaned_tfs_df, AMY_gene_summary, by = 'mgi_symbol')
  
  # Extract the number from the sheet name (e.g., 'module5' becomes '5')
  module_number <- str_extract(sheet, "\\d+")
  
  # Construct the file name using the extracted number
  output_file <- paste0("AMY_M", module_number, "_TF_DE.txt")
  
  # Write the resulting data to a file
  write.table(joinedAMY, file = output_file, row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)
  
  # Print the file name for reference
  cat("Exported:", output_file, "\n")
}