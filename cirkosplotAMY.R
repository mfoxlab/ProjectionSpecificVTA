# Set working directory
workdir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA/AMY"
colordir = "/Users/hsohail/OneDrive - Penn State Health/Documents/WGCNA"
setwd(workdir)
load(file = paste0(workdir,"/cirkosplotAMYtestworkspace.Rdata"))

# Load required libraries
library(circlize)
library(readxl)
library(RColorBrewer)
library(ComplexHeatmap)  #From BioMart

# Load the module information and module color data
module_info <- read.table("AMY_cirkosdata.txt", sep = '\t', header = TRUE, quote = '')
color_data <- read.delim(file.path(colordir, "modulecolors.txt"), header = TRUE, stringsAsFactors = FALSE, quote = '', fill = TRUE)

# Prepare the color vector by matching module names to their respective colors
module_colors <- setNames(color_data$Color, color_data$Module)

# Log base 10 transform columns 7-9 (P values)
module_info_logp <- module_info
module_info_logp[, 7:9] <- -log10(module_info[, 7:9])
colnames(module_info_logp)[7:9] <- c("logP.Value_AMY_FENT", "logP.Value_AMY_F_FENT", "logP.Value_AMY_M_FENT")
module_info_logp <- module_info_logp[order(module_info_logp$logFC_AMY_FENT), ]

# Set the specific order in which you want to plot your modules
desired_order <- module_info_logp$ModuleEigengene
# Ensure the 'modules' are ordered accordingly
modules <- factor(module_info_logp$ModuleEigengene, levels = desired_order) # WHAT DO TO HERE? DUPLICATE OF M5?
# Now subset the `module_colors` to only include the colors for `desired_order`
ordered_module_colors <- module_colors[desired_order]
# Sort module_info_logp based on desired_order to ensure alignment
module_info_logp <- module_info_logp[match(desired_order, module_info_logp$ModuleEigengene), ]


circos.clear()
pdf("cirkosplot_AMY.pdf", width = 9.5, height = 7)
circos.par(start.degree = -180)

# Initialize Circos plot with sectors in the desired order
circos.initialize(factors = modules, xlim = c(0, 1))

# Create a track and assign colors based on the specified order
circos.track(sectors = modules, ylim = c(0, 1), track.height = 0.09, 
             bg.col = module_colors[as.character(modules)],  # Apply colors in the correct order
             panel.fun = function(x, y) {
               sector_index <- CELL_META$sector.index
               circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                           labels = sector_index, 
                           facing = "inside", cex = 1,1)
             })


# Create color gradients for logFC and -log10(P value)
logFC_colors <- colorRamp2(c(-1, 0, 1), c("cyan", "white", "pink"))

# P-value gradient based on -log10(P value) scale
# Solid white for 0 to 1.3, yellow to springgreen from 1.3 to 3.5
pval_colors <- colorRamp2(c(0, 1.3, 3.5), c("white", "white", "springgreen"))

# Plot the first three tracks for logFC columns (columns 4, 5, 6)
for (i in 4:6) {
  circos.track(sectors = module_info_logp$ModuleEigengene, ylim = c(0, 1), track.height = 0.1,
               bg.col = sapply(module_info_logp[, i], logFC_colors),  # Apply logFC gradient
               panel.fun = function(x, y) {})
}

# Plot the next three tracks for P-value columns (-log10(P value) transformed in columns 7, 8, 9)
for (i in 7:9) {
  circos.track(sectors = module_info_logp$ModuleEigengene, ylim = c(0, 1), track.height = 0.1,
               bg.col = sapply(module_info_logp[, i], pval_colors),  # Apply P-value gradient
               panel.fun = function(x, y) {})
}

# Create legends for logFC and -log10(P value)
logFC_legend <- Legend(title = "logFC", col_fun = logFC_colors, 
                       at = c(-1, 0, 1), 
                       labels = c("-1", "0", "1"), 
                       title_position = "topcenter",
                       direction = "horizontal")

# Update the P-value legend to reflect the -log10(P value) scale
pval_legend <- Legend(title = "-log10(P-value)", col_fun = pval_colors, 
                      at = c(0, 1.3, 3.5),  # Corresponding to white-white-springgreen gradient
                      labels = c("0", "1.3", "3.5"), 
                      title_position = "topcenter",
                      direction = "horizontal")

# Draw legends next to the Circos plot
draw(logFC_legend, x = unit(0.8, "npc"), y = unit(0.9, "npc"), just = c("left", "top"))
draw(pval_legend, x = unit(0.8, "npc"), y = unit(0.8, "npc"), just = c("left", "top"))

dev.off()


# Clear the Circos plot space
circos.clear()


# Save the workspace
save.image(file = paste0(workdir, "/cirkosplotAMYtestworkspace.Rdata"))
