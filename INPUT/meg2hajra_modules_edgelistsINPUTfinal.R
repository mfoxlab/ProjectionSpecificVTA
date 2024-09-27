# HAJRA EDGE LIST TRIAL



# LOAD PACKAGES -----------------------------------------------------------
install.packages("reshape2")
library(reshape2)
library(dplyr)
library(WGCNA)
allowWGCNAThreads()

# PATH DIRECTORY, LOAD FILES ----------------------------------------------
path_to_files <- "C://Users/hsohail/Documents/WGCNA/INPUT" ## this is my current path but it may change for you
file_prefix <- "9TOM.INPUT" ## change this to the name of your TOM file
files <- list.files( path = path_to_files, pattern = paste0(file_prefix, "-block.*\\.RData$"), full.names = TRUE )

load(file= "C://Users/hsohail/Documents/WGCNA/INPUT/netINPUT_9.Rdata")


# FILE EDITING ------------------------------------------------------------
gene_names<-names(netINPUT_9$colors)

# Initialize an empty list to store edge lists
all_edges <- list()


for (file in files) {
  # Load the RData file
  load(file)
  
  # Convert 'dist' object to matrix
  if (exists("TOM")) {
    TOM_matrix <- as.matrix(TOM)
    diag(TOM_matrix) <- NA
    
    # Convert matrix to edge list and remove low weight
    TOM_melted <- melt(TOM_matrix, na.rm = TRUE)
    colnames(TOM_melted) <- c("Source", "Target", "Weight")
    TOM_melted <- TOM_melted[TOM_melted$Weight > 0.1, ]
    
    
    # Remove self-loops
    TOM_melted <- TOM_melted[TOM_melted$Source != TOM_melted$Target, ]
    
    # Append the edge list to the list
    all_edges[[basename(file)]] <- TOM_melted
  } else {
    warning(paste("No TOM matrix found in file:", file))
  }
}


# COMBINE EDGES -----------------------------------------------------------
# Combine all edge lists into a single data frame
combined_edges <- do.call(rbind, all_edges)
combined_edges$Source <- gene_names[as.numeric(combined_edges$Source)] 
combined_edges$Target <- gene_names[as.numeric(combined_edges$Target)]

# Add a column to the combined_edges to indicate that these are edges
combined_edges$Type <- "Edge"

# Export the combined data frame to a CSV file
write.csv(combined_edges, "C://Users/hsohail/Documents/WGCNA/INPUT/INPUT_edges.csv", row.names = FALSE)



## meg to hajra-- you can load your previous workspace, but please start here to make sure that the variables become correct 
module_assignments <- netINPUT_9$colors

node_list <- data.frame(
  ensembl_gene_id = names(module_assignments),
  ModuleID = module_assignments,
  stringsAsFactors = FALSE
)

# Extract unique gene names for node list
nodes <- unique(c(combined_edges$Source, combined_edges$Target))

# Create a data frame for the node list
nodes_df <- data.frame(Node = nodes, Type = "Node", stringsAsFactors = FALSE)

names(nodes_df) <- c("Source", "Type")

# Create an empty 'Target' and 'Weight' column to match the structure of combined_edges nodes_df$Target <- NA 
nodes_df$Weight <- NA
nodes_df$Target <- NA

#rearrange to have correct order 
nodes_df <- nodes_df[, c("Source", "Target", "Weight", "Type")]

# Combine the node list and edge list into a single data frame
combined_df <- rbind(nodes_df, combined_edges) #4.8 MILLION ROWS


## create smaller files for cytoscape to handle, based on modules

combined_edges_with_source_module <- merge(combined_df, node_list, by.x = "Source", by.y = "ensembl_gene_id") 
colnames(combined_edges_with_source_module)[ncol(combined_edges_with_source_module)] <- "Source_Module"

combined_edges_with_modules <- merge(combined_edges_with_source_module, node_list, by.x = "Target", by.y = "ensembl_gene_id")

# Rename the columns to avoid confusion
colnames(combined_edges_with_modules)[ncol(combined_edges_with_modules)] <- "Target_Module"

# Get unique modules
modules <- unique(c(combined_edges_with_modules$Source_Module, combined_edges_with_modules$Target_Module))
# Directory to save files
output_dir <- "module_edges_input"
dir.create(output_dir, showWarnings = FALSE)
# Loop through each module and save the corresponding edges
for (module in modules) {
  # Filter edges where either source or target belongs to the module
  module_edges <- subset(combined_edges_with_modules, Source_Module == module & Target_Module == module)
  
  # File name based on the module
  file_name <- paste0(output_dir, "/module_", module, "_edges.txt")
  
  # Save the edges
  write.table(module_edges, file_name, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Print a message
  cat("Saved edges for module", module, "to", file_name, "\n")
}

save.image(file = "C://Users/hsohail/Documents/WGCNA/INPUT/INPUTEdgesWorkspace.RData")


