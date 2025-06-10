# Set up differential expression and remove genes that do not have at least 1 CPM in at least 2 samples   

# SETUP ------------------------------------------------------------------------------------------

# 1: Download counts.Rdata and meta.Rdata from Github
# 2: Set your working directory in line 7

workdir = ""
setwd(workdir)

library(limma)
library(edgeR)
library(data.table)

load("counts.Rdata")
load("meta.Rdata")


# QUALITY CONTROL ------------------------------------------------------------------------------------------

y <- DGEList(counts = counts, group = meta$Group) 

keep <- rowSums(cpm(y) > 1) >= 2 & grepl("ENSMUSG", rownames(cpm(y))) 

y = y[keep, , keep.lib.sizes = T] 
y <- calcNormFactors(y) 
logcpm = cpm(y, log = T, prior.count = 3)


#make MDS plots 
colnames(y$counts)<-y$samples$group
col = vector(mode = "character", length = length(colnames(y$counts)))
pch = vector(mode = "character", length = length(colnames(y$counts)))

pch[grep("FF", colnames(y$counts))] = 0
pch[grep("MF", colnames(y$counts))] = 16
pch[grep("FS", colnames(y$counts))] = 15
pch[grep("MS", colnames(y$counts))] = 1


col[grep("AMY", colnames(y$counts))] = "hotpink"
col[grep("NAC", colnames(y$counts))] = "red2"
col[grep("INPUT", colnames(y$counts))] = "blue"
col[grepl("PFC", colnames(y$counts))] = "springgreen3"

#everything combined
pdf("mds_allsamples.pdf")
plotMDS(y, dim.plot=c(1,2), cex=1.2, top = 500, col = col, pch= as.numeric(pch))
legend("topright", col = c("springgreen3", "red2","blue"," hotpink"), legend=c("PFC", "NAc", "Input", "AMY"), pch=20, cex=1)
legend("bottomright", col = c("grey77"), legend=c("MFent", "MSal","FFent","FSal"), pch= c(16, 1,0, 15), cex = 1)
dev.off()

#just VTA->PFC
pdf("mds_pfc.pdf")
plotMDS(y[, meta$target == "PFC"], dim.plot=c(1,2),col="springgreen3", cex=1.2, pch=as.numeric(pch[meta$target == "PFC"]), main = "PFC")
legend("bottomleft", col = c("springgreen3"), legend=c("MFent", "MSal", "FFent","FSal"),pch= c(16,1,0,15),cex=1)
dev.off()


#just VTA->AMY
pdf("mds_amy.pdf")
plotMDS(y[, meta$target == "AMY"], dim.plot=c(1,2),col="hotpink", cex=1.2, pch=as.numeric(pch[meta$target == "PFC"]), main = "PFC")
legend("bottomleft", col = c("hotpink"), legend=c("MFent", "MSal", "FFent","FSal"),pch= c(16,1,0,15),cex=1)
dev.off()

#just VTA->NAC     
pdf("mds_nac.pdf")
plotMDS(y[, meta$target == "NAC"], dim.plot=c(1,2),col="red2", cex=1.2, pch=as.numeric(pch[meta$target == "NAC"]), main = "NAc")
legend("bottomleft", col = c("red2"), legend=c("MFent", "MSal", "FFent","FSal"),pch= c(16,1,0,15),cex=1)
dev.off()


#just VTA “input” samples  
pdf("mds_input.pdf")
plotMDS(y[, meta$target == "INPUT"], dim.plot=c(1,2),col="blue", cex=1.2, pch=as.numeric(pch[meta$target == "INPUT"]), main = "INPUT")
legend("bottomleft", col = c("blue"), legend=c("MFent", "MSal", "FFent","FSal"),pch= c(16,1,0,15),cex=1)
dev.off()

# DIFFERENTIAL EXPRESSION WITH lmfit ------------------------------------------------------------------------------------------
unique_groups <- unique(y$samples$group)

design <- matrix(0, nrow = nrow(y$samples), ncol = length(unique_groups))

design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
for (i in 1:nrow(y$samples)) {
  design[i, y$samples$group[i]] <- 1
}

fit1 = lmFit(logcpm, design = design) 
fit1 <- eBayes(fit1)

#comparisons within brain regions, separately in male and female, fent vs saline, as well as sexes combined, plus brain regions vs input, plus brain regions vs each other

contr = makeContrasts(
  NAC_F_FENT=	(NACFFENT-NACFSAL),
  NAC_M_FENT=	(NACMFENT-NACMSAL),
  NAC_FENT=	(NACFFENT+NACMFENT)/2-(NACFSAL+NACMSAL)/2,
  NAC_FENT_INPUT=	(NACFFENT+NACMFENT)/2 - (INPUTFFENT+INPUTMFENT)/2,
  NAC_SAL_INPUT=	(NACFSAL+NACMSAL)/2 - (INPUTFSAL+INPUTMSAL)/2,
  PFC_F_FENT=	(PFCFFENT-PFCFSAL),
  PFC_M_FENT=	(PFCMFENT-PFCMSAL),
  PFC_FENT=	(PFCFFENT+PFCMFENT)/2-(PFCFSAL+PFCMSAL)/2,
  PFC_FENT_INPUT=	(PFCFFENT+PFCMFENT)/2-(INPUTFFENT+INPUTMFENT)/2,
  PFC_SAL_INPUT=	(PFCFSAL+PFCMSAL)/2-(INPUTFSAL+INPUTMSAL)/2,
  INPUT_F_FENT=	(INPUTFFENT-INPUTFSAL),
  INPUT_M_FENT=	(INPUTMFENT-INPUTMSAL),
  INPUT_FENT=	(INPUTFFENT+INPUTMFENT)/2-(INPUTFSAL+INPUTMSAL)/2, 
  NAC_PFC_FENT = (NACFFENT+NACMFENT)/2-(PFCFFENT+PFCMFENT)/2,
  NAC_PFC_SAL = (NACFSAL+NACMSAL)/2-(PFCFSAL+PFCMSAL)/2,
  levels=design)

fit1a = contrasts.fit(fit1, contr) 
fit1a = eBayes(fit1a, trend = T) 
res1 = topTable(fit1a, coef = 1, number = Inf) 
res1 = merge(res1, topTable(fit1a, coef = 2, number = Inf), by = 0, suffixes = c(".coef1", ".coef2")) 
res1 = merge(res1, topTable(fit1a, coef = 3, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef3")) 
res1 = merge(res1, topTable(fit1a, coef = 4, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef4")) 
res1 = merge(res1, topTable(fit1a, coef = 5, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef5"))
res1 = merge(res1, topTable(fit1a, coef = 6, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef6"))
res1 = merge(res1, topTable(fit1a, coef = 7, number = Inf), by.x ="Row.names", by.y = 0, suffixes = c("", ".coef7"))
res1 = merge(res1, topTable(fit1a, coef = 8, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef8"))
res1 = merge(res1, topTable(fit1a, coef = 9, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef9"))
res1 = merge(res1, topTable(fit1a, coef = 10, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef10"))
res1 = merge(res1, topTable(fit1a, coef = 11, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef11"))
res1 = merge(res1, topTable(fit1a, coef = 12, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef12"))
res1 = merge(res1, topTable(fit1a, coef = 13, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef13"))
res1 = merge(res1, topTable(fit1a, coef = 14, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef14"))
res1 = merge(res1, topTable(fit1a, coef = 15, number = Inf), by.x = "Row.names", by.y = 0, suffixes = c("", ".coef15"))

n = ncol(contr) 
limma.pvals = res1[, c(1:n * 6 - 1)]   #extract p value
limma.logFC = res1[, c(1:n * 6 - 4)]   #extract log fold change
limma.AveExpr = res1[, c(1:n * 6 - 3)] #extract average expression
colnames(limma.pvals) = colnames(limma.logFC) = colnames(limma.AveExpr) = colnames(contr) 
rownames(limma.pvals) = rownames(limma.logFC) = rownames(limma.AveExpr) = res1[, 1] 

limma.fdr = apply(limma.pvals, 2, p.adjust) 
colnames(limma.fdr) = paste("FDR", colnames(limma.fdr), sep = "_") 
colnames(limma.pvals) = paste("P.Value", colnames(limma.pvals), sep = "_") 
colnames(limma.logFC) = paste("logFC", colnames(limma.logFC), sep = "_") 
colnames(limma.AveExpr) = paste("AveExpr", colnames(limma.AveExpr), sep = "_") 

summa.fit1a <- decideTests(fit1a, adjust.method = "none", p.value = 0.001)
summary_table <- data.frame(Contrast = character(), Down = numeric(), NotSig = numeric(), Up = numeric(), stringsAsFactors = FALSE)

for (col in colnames(summa.fit1a)) {
  #extract decision results for the current contrast
  decisions <- summa.fit1a[, col]
  
  down_count <- sum(decisions == -1)
  notsig_count <- sum(decisions == 0)
  up_count <- sum(decisions == 1)
  
  summary_table <- rbind(summary_table, data.frame(Contrast = col, Down = down_count, NotSig = notsig_count, Up = up_count))
}

#annotate with mgi_symbols
library(biomaRt)
mart = useMart("ensembl") 
mart = useDataset("mmusculus_gene_ensembl", mart) 
ann = getBM(mart = mart, attributes = c("ensembl_gene_id", "mgi_symbol", "gene_biotype", "description"), filters = "ensembl_gene_id", values = rownames(logcpm))
ann = ann[duplicated(ann$ensembl_gene_id) == F, ]

target=unique(meta$target)

#based on p<0.05 
sig.genes.pval=list()
pcut=0.05
for (i in target) {sig.genes.pval[[i]] = rownames(limma.pvals)[rowSums(limma.pvals[,grep(i,colnames(limma.pvals))] <= pcut)>0] }

#based on FDR <0.1
sig.genes.fdr=list()
for(i in target) { sig.genes.fdr[[i]] = rownames(limma.fdr)[rowSums(limma.fdr[,grep(i,colnames(limma.fdr))] <=0.1)>0]}
cbind(genepval=lapply(sig.genes.pval,length),genefdr=lapply(sig.genes.fdr,length))

'''
results here
      genepval genefdr
PFC   18801    9658   
AMY   0        0      
INPUT 18888    11169  
NAC   18666    10352  
'''

#write  genes p<0.05
sigGenes.tab = lapply(1:length(target), function(i) { data.frame(ensembl_gene_id = sig.genes.pval[[target[i]]], limma.pvals[sig.genes.pval[[target[i]]], grep(target[i], colnames(limma.pvals))], limma.logFC[sig.genes.pval[[target[i]]], grep(target[i], colnames(limma.pvals))]) }) 
names(sigGenes.tab) = target 
DEoutp = lapply(1:length(target), function(i) { merge(ann, sigGenes.tab[[target[i]]], by = 1) }) 
names(DEoutp) = target 

save(DEoutp, file = "DEoutputp05.Rdata")

#create individual lists of DEGs per grouping
for (i in 1:length(target)) {
  flname <- paste0(getwd(), "/DEgene_summary_", target[i], ".txt")
  write.table(DEoutp[[target[i]]], file = flname, sep = "\t", row.names=FALSE)
}
library(dplyr) 
#create the bio_rep column based on the incidence within Group for easier organization later
meta$bio_rep <- ave(seq_along(meta$Group), meta$Group, FUN = seq_along)

#create new sample names using target, drug, sex, and bio_rep
new_sample_names <- paste(meta$target, meta$drug, meta$sex, meta$bio_rep, sep = "_")

#rename columns in the original logcpm matrix
colnames(logcpm) <- new_sample_names

#subset logcpm for each region based on significant genes
logcpm_PFC_sig <- logcpm[rownames(logcpm) %in% sig.genes.pval$PFC, grep("PFC", colnames(logcpm))]
logcpm_NAC_sig <- logcpm[rownames(logcpm) %in% sig.genes.pval$NAC, grep("NAC", colnames(logcpm))]
logcpm_INPUT_sig <- logcpm[rownames(logcpm) %in% sig.genes.pval$INPUT, grep("INPUT", colnames(logcpm))]


#save the logCPM objects to use for further analysis
save(logcpm_PFC_sig, file = "logcpm_PFC_sig.Rdata")
save(logcpm_NAC_sig, file = "logcpm_NAC_sig.Rdata")
save(logcpm_INPUT_sig, file = "logcpm_INPUT_sig.Rdata")


