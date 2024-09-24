
## load your file here (i did this)
workdir = "/Users/hsohail/Documents/WGCNA"
workdirPFC = "/Users/hsohail/Documents/WGCNA/PFC"
setwd(workdirPFC)

load(paste0(workdir, "/logcpm_PFC_sig.Rdata"))

## building WGCNA using the data with only DEGs as function of drug or projection target. Example shown for VTA->PFC

library(WGCNA)
allowWGCNAThreads()


# Define the power range 
powers <- c(1:30, seq(from = 35, to = 50, by = 5)) 

#transpose the data so that genes and samples are correctly arrangedd
#logcpm_PFC_sig<- t(logcpm_PFC_sig)

# Perform the soft-thresholding analysis for VTA->PFC
powerTable <- pickSoftThreshold(logcpm_PFC_sig, powerVector = powers, verbose = 5)

# Plot the results for this set
par(mfrow = c(1, 2)) 
plot(powerTable$fitIndices[, 1], -sign(powerTable$fitIndices[, 3]) * powerTable$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
     main = paste("Scale Independence for PFC"))
text(powerTable$fitIndices[, 1], -sign(powerTable$fitIndices[, 3]) * powerTable$fitIndices[, 2], labels = powers, col = "red") 
plot(powerTable$fitIndices[, 1], powerTable$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity for PFC")) 
text(powerTable$fitIndices[, 1], powerTable$fitIndices[, 5], labels = powers, col = "red")


#find the modules using blockwise

powerthres = 9
netPFC_9 = blockwiseModules(logcpm_PFC_sig, corType = "bicor", maxPOutliers = 0.1, 
                            power = powerthres, networkType = "signed", 
                            minModuleSize = 30, reassignThreshold = 0,
                            mergeCutHeight = 0.25, numericLabels = TRUE, minMEtoStay = 0, pamRespectsDendro=FALSE, 
                            saveTOMs = TRUE, saveTOMFileBase = "9TOM.PFC",
                            loadTOM = T, verbose = 5)

save(netPFC_9,file="netPFC_9.Rdata")

#relate modules to traits
sample_names<-rownames(logcpm_PFC_sig)
traits <-ifelse(grepl("FENT", sample_names), 1, 0)
trait_data<-data.frame(Sample=sample_names, Trait=traits)
rownames(trait_data)<-sample_names

#define module colors
PFC_module_colors<- netPFC_9$colors

# look at expression as function of fentanyl treatment
design <- model.matrix(~ trait_data$Trait)
expr_data<- t(logcpm_PFC_sig)
results <- list()
module_list <- unique(PFC_module_colors) 


results <- list()
library(limma)
for(module in module_list) { 
  module_genes <- names(PFC_module_colors)[PFC_module_colors == module] 
  module_expr_data <- expr_data[rownames(expr_data) %in% module_genes, ] 
  gene_index <- which(rownames(expr_data) %in% module_genes) 
  roast_result <- roast(expr_data, index = gene_index, design = design)
  results[[as.character(module)]] <- roast_result
}

#identify MEs and make adjacency heatmaps

#identify module eigengenes for PFC

MEsPFC=moduleEigengenes(logcpm_PFC_sig, PFC_module_colors, impute=TRUE, nPC=1, align= "along average", excludeGrey=FALSE, grey = if( is.numeric(PFC_module_colors))0 else "grey", subHubs=TRUE,softPower=9, scale=TRUE, verbose=5, indent=0)

save(MEsPFC, file= "MEsPFC.RData")

library(pheatmap) 
eigengenes <- MEsPFC$eigengenes
adjacency_matrix <- cor(eigengenes, use = "pairwise.complete.obs")


# Plot the adjacency matrix using pheatmap
pdf("pfc_adjacency_matrix_heatmap.pdf")
pheatmap(adjacency_matrix,
         main = "Adjacency Matrix of Module Eigengenes_PFC",
         color = colorRampPalette(c("blue", "white", "red"))(50), # Customize the color palette
         cluster_rows = TRUE, # Cluster the rows
         cluster_cols = TRUE, # Cluster the columns
         display_numbers = FALSE)

dev.off()

# save/ export the file with module membership for all the genes

# load in the DE tables to append to the gene summary 
load("DEgene_summary_PFC.txt") #HAJRA NEW CHECK THIS

DEsPFC<- read.table("DEgene_summary_PFC.txt", sep= "\t", header=TRUE)


geneModuleMembership = list()
MEsPFCAve = MEsPFC$averageExpr
nSamples= nrow(logcpm_PFC_sig)
ModuleMembership= as.data.frame(bicor(logcpm_PFC_sig, MEsPFCAve, use = "p"))
MMPValue = as.data.frame(corPvalueStudent(as.matrix(ModuleMembership), nSamples))
names(ModuleMembership) = gsub("AE", paste0("PFC", ".MM"), names(ModuleMembership))
names(MMPValue) = gsub("AE", "p.MM", names(ModuleMembership))
geneModuleMembership = list(membership=as.matrix(ModuleMembership), pvals=as.matrix(MMPValue))

save(geneModuleMembership, file="PFC_gene_module_membership.Rdata") 

netGenes= data.frame(ensembl_gene_id = colnames(logcpm_PFC_sig), module.PFC = netPFC_9$colors)


# Convert membership matrix to data frame for easy joining 
membership_df <- as.data.frame(geneModuleMembership$membership, stringsAsFactors = FALSE)
membership_df$ensembl_gene_id <- rownames(membership_df)

library(dplyr)
outp <- DEsPFC %>%
  full_join(netGenes, by = "ensembl_gene_id") %>%
  full_join(membership_df, by = "ensembl_gene_id") 


# Export the data 
write.table(outp, "PFC_gene_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# find hub genes 
library(biomaRt) #ADDED AFTER RUNNING
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

topHubGenesList <- list()

moduleNames <- colnames(ModuleMembership)

for (module in moduleNames) {
  moduleMembershipValues <- ModuleMembership[[module]]
  
  rankedGenes <- order(moduleMembershipValues, decreasing = TRUE)
  
  
  top10Genes <- rankedGenes[1:10]
  
  
  top10GeneIDs <- rownames(ModuleMembership)[top10Genes]
  top10kME <- moduleMembershipValues[top10Genes]
  
  geneModuleMembership<-as.data.frame(geneModuleMembership)
  annotations <- getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(geneModuleMembership),
    mart = ensembl
  )
  
  
  topHubGenesDF <- data.frame(
    ensembl_gene_id = top10GeneIDs,
    kME = top10kME
  )
  topHubGenesDF <- merge(topHubGenesDF, annotations, by = "ensembl_gene_id", all.x = TRUE)
  
  topHubGenesList[[module]] <- topHubGenesDF
}
write.table(topHubGenesList, "tophubgeneslist_pfc.txt", sep="\t", row.names=FALSE)




###
  #alternative strategy to determine fentanyl associated modules 
  
  
  combined_df<-cbind(MEsPFC$eigengenes, Condition=trait_data$Trait)
combined_df$Condition <- factor(combined_df$Condition, levels = c(0, 1), labels = c('SAL', 'FENT'))
kME <- signedKME(logcpm_PFC_sig, MEsPFC$eigengenes)
hub_genes <- apply(kME, 2, function(x) names(x)[which.max(x)])
names(hub_genes) <- paste0("ME", sub("kME", "", names(hub_genes)))


# Create design matrix

design <- model.matrix(~ Condition, data = combined_df)
eigengenes_df <- t(MEsPFC$eigengenes)

fit <- lmFit(eigengenes_df, design)
fit <- eBayes(fit)

# Get the results table
results <- topTable(fit, coef="ConditionFENT", number=nrow(eigengenes_df)) 
results$ModuleEigengene <- rownames(results)
results$ensembl_gene_id <- sapply(rownames(results), function(x) hub_genes[x]) 

#annotate with mgi_symbols
results_annotated <- merge(results, annotations, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x=TRUE)

##note ann is from a previous workspace-- used because of biomart problemes #HAJRA EDITED UNCOMMENTED 192-196 DO THESE LINES FIRST BEFORE ROW 189



write.table(results_annotated, "PFCmoduleEigengenesDE.txt", sep="\t", row.names=FALSE, quote = FALSE)

save.image(file = "logCPM_PFC_workspace.RData")
