#script to run DESeq2
#10/20/2020 MF

##Install DESeq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("DESeq2")
#BiocManager::install('EnhancedVolcano')

#Load Libraries
library(DESeq2); library(pheatmap); library(ashr); library(EnhancedVolcano); library(magrittr); library(ggfortify); library(reshape2); library(ggplot2); library(GOplot); library(arm)

#set working directory - tried multiple quality scores and aligning by gene and exon (no real difference in output)
#so sticking with hiQual_ex
hiQual_Ex<-("/media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_exon/")
#loQual_Ex<-("/media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/lowQual_exon/")
#hiQual_Ge<-("/media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_gene/")
#loQual_Ge<-("/media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/lowQual_gene/")

setwd(hiQual_Ex)

#Load in read counts and assign them sample labels
outputPrefix<-("Culex_DEseq")
sampleFiles<-c("M1-1_S1_htseq","M1-2_S2_htseq","M1-3_S3_htseq","M1-4_S4_htseq","M2-1_S5_htseq","M2-2_S6_htseq","M2-3_S7_htseq","M2-4_S8_htseq","M4-1_S13_htseq","M4-2_S14_htseq","M4-3_S15_htseq","M4-4_S16_htseq")
sampleNames<-c("BG_Gravid1", "BG_Gravid2", "BG_Gravid3", "BG_Gravid4", "BG_Parous1", "BG_Parous2", 
               "BG_Parous3", "BG_Parous4", "AG1", "AG2", "AG3", "AG4")
sampleCondition<-c("CALgravidF",        "CALgravidF",   "CALgravidF",   "CALgravidF",   "CALparousF",   "CALparousF",   "CALparousF",   "CALparousF",  "PipEvanF",     "PipEvanF",     "PipEvanF",     "PipEvanF")
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)
treatments<-c("CALgravidF",     "CALparousF",   "PipEvanF")

#Create DESeq Data - CHANGE DIRECTORY
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = hiQual_Ex,
                                       design= ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels=treatments)
dim(ddsHTSeq)#getting num genes in dataset and verifying num samples


#Prefilter so only genes with at least 10 reads in at least 4 samples are considered
keep <- rowSums(counts(ddsHTSeq) >= 10) >= 4
ddsHTSeq <- ddsHTSeq[keep,]
dim(ddsHTSeq)#new filtered num genes in dataset and verifying num samples

ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = treatments)



#Looking at distr of filtered count data across samples
librarySizes <- colSums(counts(ddsHTSeq))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(ddsHTSeq) + 1)
head(logcounts)

#Any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(ddsHTSeq$condition)) + 1  # make a colour vector

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")



#Looking at PCA of the data - Do treatments cluster together?
rlogcounts <- rlog(counts(ddsHTSeq))#transforming data to make it approximately homoskedastic, n < 30 so rlog is better

select = order(rowMeans(rlogcounts), decreasing=TRUE)[1:12710]
#This select variable is here because I modified the numbers of genes included in the PCA - from 500-12710.
#The percent variation explained by each PC changes with the number of genes included, but not THAT much.
#At lower numbers of genes, PhysStat can be predicted by PC1 better than PC2, but this changes as I increase to include more.
#Because I couldn't decide on a sensible cutoff for the number of genes to include, I used the full dataset.


highexprgenes_counts <- rlogcounts[select,]

colnames(highexprgenes_counts)<- ddsHTSeq$condition

data_for_PCA <- t(highexprgenes_counts)
dim(data_for_PCA)

#run PCA
pcDat <- prcomp(data_for_PCA, center = T)

#basic plot
autoplot(pcDat)

#plot for pub
pdf("Fig2_PCA_Treatments.pdf",width=6,height=6,paper='special')
autoplot(pcDat,
         data = ddsHTSeq$colData, 
         colour=as.numeric(factor(ddsHTSeq$condition)),
         shape=FALSE, 
         label.size=6, xlim = c(-0.4, 0.5)) + theme_bw()
dev.off()


#test of whether PCs can predict strain and physState
Treats <- as.character(c("BG_grav", "BG_grav", "BG_grav", "BG_grav", "BG_par", "BG_par", "BG_par", "BG_par", "AG2", "AG2", "AG2", "AG2"))
PC1 <- as.character(pcDat$x[,1])
PC2 <- as.character(pcDat$x[,2])
PC3 <- as.character(pcDat$x[,3])
Strain <- as.character(c("BG", "BG", "BG", "BG", "BG", "BG", "BG", "BG", "AG", "AG", "AG", "AG"))
PhysStat <- as.character(c("grav", "grav", "grav", "grav", "hs", "hs", "hs", "hs", "hs", "hs", "hs", "hs"))

dat_for_glm <- data.frame(cbind(Treats, Strain, PhysStat, PC1, PC2, PC3), row.names = NULL)
summary(dat_for_glm)
str(dat_for_glm)

dat_for_glm$PC1 <- as.numeric(as.character(dat_for_glm$PC1))
dat_for_glm$PC2 <- as.numeric(as.character(dat_for_glm$PC2))
dat_for_glm$PC3 <- as.numeric(as.character(dat_for_glm$PC3))


#Are PCs 1,2,3 capable of separating samples by strain?
Model_strain_PC1 <- bayesglm(Strain ~ PC1, data=dat_for_glm, family="binomial")
summary(Model_strain_PC1)

simulates <- coef(sim(Model_strain_PC1, n.sims = 10000))
head(simulates, 10)
plot(density(simulates[,2]), main = "Model_strain_PC1", xlab = "Posterior.open", ylab = "Density")

quantile(simulates[,2], c(0.025, 0.975))#gives the 95% credible intervals, which don't overlap with zero
#This indicates that PC1 is capable of separating our samples by strain.

Model_strain_PC2 <- bayesglm(Strain ~ PC2, data=dat_for_glm, family="binomial")
summary(Model_strain_PC2)

simulates <- coef(sim(Model_strain_PC2, n.sims = 10000))
head(simulates, 10)
plot(density(simulates[,2]), main = "Model_strain_PC2", xlab = "Posterior.open", ylab = "Density")

quantile(simulates[,2], c(0.025, 0.975))#PC2 - overlaps zero

Model_strain_PC3 <- bayesglm(Strain ~ PC3, data=dat_for_glm, family="binomial")
summary(Model_strain_PC3)

simulates <- coef(sim(Model_strain_PC3, n.sims = 10000))
head(simulates, 10)
plot(density(simulates[,2]), main = "Model_strain_PC3", xlab = "Posterior.open", ylab = "Density")

quantile(simulates[,2], c(0.025, 0.975))#PC3 - overlaps zero


#Are PCs 1,2,3 capable of separating samples by physState?
Model_PhysStat_PC1 <- bayesglm(PhysStat ~ PC1, data=dat_for_glm, family="binomial")
summary(Model_PhysStat_PC1)

simulates <- coef(sim(Model_PhysStat_PC1, n.sims = 10000))
head(simulates, 10)
plot(density(simulates[,2]), main = "Model_PhysStat_PC1", xlab = "Posterior.open", ylab = "Density")

quantile(simulates[,2], c(0.025, 0.975))#PC1 & PhysStat - overlaps zero.

Model_PhysStat_PC2 <- bayesglm(PhysStat ~ PC2, data=dat_for_glm, family="binomial")
summary(Model_PhysStat_PC2)

simulates <- coef(sim(Model_PhysStat_PC2, n.sims = 10000))
head(simulates, 10)
plot(density(simulates[,2]), main = "Model_PhysStat_PC2", xlab = "Posterior.open", ylab = "Density")

quantile(simulates[,2], c(0.025, 0.975))#PC2 & PhysStat - does not overlap zero

Model_PhysStat_PC3 <- bayesglm(PhysStat ~ PC3, data=dat_for_glm, family="binomial")
summary(Model_PhysStat_PC3)

simulates <- coef(sim(Model_PhysStat_PC3, n.sims = 10000))
head(simulates, 10)
plot(density(simulates[,2]), main = "Model_PhysStat_PC3", xlab = "Posterior.open", ylab = "Density")

quantile(simulates[,2], c(0.025, 0.975))#PC3 & PhysStat - overlaps zero



#Note that the PCA shows that BG_Parous4 is fairly different from the other BG parous treatments - should I drop it?  See end of script for code to do this - in the end, I did not for the paper.


####Differential Expression Analysis
##First with all data
dds<-DESeq(ddsHTSeq)
res<-results(dds)

resultsNames(dds)

###getting normalized filtered read counts
dds <- estimateSizeFactors(dds)
Norm_counts <- counts(dds, normalized=TRUE)
cormat <- cor(Norm_counts, method = "spearman")#looking at correlation coefficients for gene expression values between treatments.

#Function to Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "pink", 
                       midpoint = 0.95, limit = c(0.90,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  xlab("") + ylab("")+
  coord_fixed()

# Print the heatmap
print(ggheatmap)

#write.csv(as.data.frame(Norm_counts), 
          #file="Normalized_ReadCounts_AllTreats.csv")
#Norm_Rownames <- data.frame(rownames(Norm_counts))
#Norm_Rownames$Descriptor <- "NA"

#CALsOnly <- cbind(Norm_Rownames, (data.frame(Norm_counts[,c(1:8)])))
#write.table(as.data.frame(CALsOnly), 
      #file="Normalized_ReadCounts_CALsOnly.txt", row.names = F, sep = "\t")

#CalParVPip<- cbind(Norm_Rownames, (data.frame(Norm_counts[,c(-1,-2,-3,-4)])))
#write.table(as.data.frame(CalParVPip), 
      #file="Normalized_ReadCounts_CalParVPip.txt", row.names = F, sep = "\t")


###Gravid F vs. Parous F--Use LFC for gene ranking and visualization, and use the p-values from the non-LFC
resLFC <- lfcShrink(dds, contrast = c("condition", "CALgravidF", "CALparousF"), type="ashr")

#Look at summary values
summary(resLFC)

#how many significantly DE genes? Looked at adjusted pvals of 0.1 and 0.05. Added a FC cutoff for shrunken FCs of 1.5x, as well.
sum(resLFC$padj < 0.1 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)
sum(resLFC$padj < 0.05, na.rm=TRUE)

sum(resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)

#getting numbers of up and down regulated genes
upReg <- subset(resLFC, padj < 0.05 & log2FoldChange > 0.58)#LFC > 1.5 or upReg in gravid
nrow(upReg)
head(upReg)

downReg <- subset(resLFC, padj < 0.05 & log2FoldChange < -0.58)#LFC < 1.5 or downReg in gravid
nrow(downReg)
head(downReg)


#Create plots based on LFCshrunken dataset, which minimizes noise from low read counts
pdf("plotMA_gravidVparous.pdf",width=6,height=6,paper='special')
plotMA(resLFC, ylim=c(-3,3))
dev.off()

#plotting differences in gene expression using lfcShrink output.
pdf("EV_gravidVparous.pdf",width=8,height=6,paper='special') #gene names after selectLab are all diff exp genes
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c('CPIJ003142', 'CPIJ004468', 'CPIJ004690', 'CPIJ008018', 'CPIJ008747', 'CPIJ010041', 'CPIJ011084', 'CPIJ011244',
                #'CPIJ015908','CPIJ018848', 'CPIJ003456', 'CPIJ004365', 'CPIJ004417', 'CPIJ008256', 'CPIJ012990', 'CPIJ014981'),
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)
dev.off()


#Looking at genes using the standard frequentist framework with a false discovery rate correction. 
res05 <- results(dds, alpha=0.05, contrast = c("condition", "CALgravidF", "CALparousF"))#results here are largely consistent with the shrunken LFC framework above without the LFC threshold argument. 
summary(res05) #note there are quite a few with "low" counts (mean count < 88), total significant = 556 compared with the 544 above.

#non-LFC dataset filtering
sum(res05$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)

# Note the small difference between LFC shrunken and "standard" analyses are normal - inconsistencies between them are typically minimal and have been described: https://support.bioconductor.org/p/110307/


#getting full gene lists with different LFCs and padj - from lfcShrink and results.  Note that genes themselves are the same, what differ (mildly) are the statistics.
full_genes_CalgravidF_CalparousF <- data.frame(c(paste(resLFC@rownames, sep = "", "-RA")), resLFC@listData$log2FoldChange, resLFC@listData$padj)#paste adds VB format to IDs
write.table(full_genes_CalgravidF_CalparousF, file = "CALgravidF_v_CALparousF_AllGenes_LFCS.txt", 
            col.names = F, row.names = F, sep = "\t")

full_genes_CalgravidF_CalparousF <- data.frame(c(paste(res05@rownames, sep = "", "-RA")), res05@listData$log2FoldChange, res05@listData$padj)#paste adds VB format to IDs
write.table(full_genes_CalgravidF_CalparousF, file = "CALgravidF_v_CALparousF_AllGenes_nonLFCS.txt", 
            col.names = F, row.names = F, sep = "\t")

#Saving a list of significant DE genes
resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) > 0.58)
write.table(data.frame(c(paste(resSig@rownames, sep = "", "-RA")), resSig@listData$log2FoldChange, resSig@listData$padj), 
            file="CALgravidF_v_CALparousF_LFCS_padj05.txt", col.names = F, row.names = F, sep = "\t")

resSig <- subset(res05, padj < 0.05 & abs(log2FoldChange) > 0.58)
write.table(data.frame(c(paste(resSig@rownames, sep = "", "-RA")), resSig@listData$log2FoldChange, resSig@listData$padj),
          file="CALgravidF_v_CALparousF_nonLFCS_padj05.txt",  col.names = F, row.names = F, sep = "\t")




###Parous F vs. pipiens
#######Here looking at full dataset including the BG_parous4 outlier

resLFC <- lfcShrink(dds, contrast = c("condition", "PipEvanF", "CALparousF"), type="ashr")

#Look at summary values
summary(resLFC)

#how many significantly DE genes? The default p-value cutoff is 0.1 & adding a log2FoldChange cutoff (1.5), as well.
sum(resLFC$padj < 0.1 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)
sum(resLFC$padj < 0.05, na.rm=TRUE)

sum(resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)

upReg <- subset(resLFC, padj < 0.05 & log2FoldChange > 0.58)#LFC > 1.5 or upReg in Pip
nrow(upReg)
head(upReg)

downReg <- subset(resLFC, padj < 0.05 & log2FoldChange < -0.58)#LFC < 1.5 or downReg in Pip
nrow(downReg)
head(downReg)

#Create plots based on the LFC, which minimizes noise from low read counts
pdf("plotMA_parousVpip.pdf",width=6,height=6,paper='special')
plotMA(resLFC, ylim=c(-3,3))
dev.off()

pdf("EV_parousVpip.pdf",width=8,height=6,paper='special') #gene labels commented out for selectLab are sig DE sensory genes.
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c("CPIJ001730", "CPIJ002108", "CPIJ002109", "CPIJ002111", "CPIJ004145","CPIJ007617",
                              #"CPIJ009568", "CPIJ010367", "CPIJ010787", "CPIJ012716","CPIJ012717", "CPIJ012719",
                             #"CPIJ013976", "CPIJ014525", "CPIJ016479","CPIJ016949", "CPIJ016966", "CPIJ019610",
                             #"CPIJ016433", "CPIJ011564", "CPIJ014330","CPIJ002605", "CPIJ002618", "CPIJ002628",
                              #"CPIJ007315", "CPIJ004067"),
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-10, 10),
                ylim = c(0,200),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 2.0)
dev.off()

#Looking at DGE from "standard" analysis.
res05 <- results(dds, alpha=0.05, contrast = c("condition", "PipEvanF", "CALparousF"))
summary(res05)
sum(res05$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)


#getting full gene list
full_genes_CalparousF_PipF <- data.frame(c(paste(resLFC@rownames, sep = "", "-RA")), resLFC@listData$log2FoldChange, resLFC@listData$padj)#paste adds VB format to IDs
write.table(full_genes_CalparousF_PipF, file = "CALparousF_v_PipEvanF_AllGenes_LFCS.txt", 
            col.names = F, row.names = F, sep = "\t")

full_genes__CalparousF_PipF <- data.frame(c(paste(res05@rownames, sep = "", "-RA")), res05@listData$log2FoldChange, res05@listData$padj)#paste adds VB format to IDs
write.table(full_genes_CalparousF_PipF, file = "CALparousF_v_PipEvanF_AllGenes_nonLFCS.txt", 
            col.names = F, row.names = F, sep = "\t")



#Save a list of significant DE genes with 1.5x change or greater
resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) > 0.58)
write.table(data.frame(c(paste(resSig@rownames, sep = "", "-RA")), resSig@listData$log2FoldChange, resSig@listData$padj),
            file="CALparousF_v_PipEvanF_LFCS_padj05.txt", row.names = F, col.names = F, sep = "\t") 

resSig <- subset(res05, padj < 0.05 & abs(log2FoldChange) > 0.58 )
write.table(data.frame(c(paste(resSig@rownames, sep = "", "-RA")), resSig@listData$log2FoldChange, resSig@listData$padj),
                        file="CALparousF_v_PipEvanF_nonLFCS_padj05.txt", row.names = F, col.names = F, sep = "\t") 
  


#####This contrast is not that interesting, but looked at it anyway. Will be useful to verify those 16 genes that look important to suppressing host-seeking in BG gravid.
###Gravid F vs. pipiens

resLFC <- lfcShrink(dds, contrast = c("condition",  "PipEvanF", "CALgravidF"), type="ashr")

#Look at summary values
summary(resLFC)

#how many significantly DE genes? The default p-value cutoff is 0.1. Added a log2foldchange cutoff of 1.5x, as well.
sum(resLFC$padj < 0.1 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)
sum(resLFC$padj < 0.05, na.rm=TRUE)

sum(resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)

upReg <- subset(resLFC, padj < 0.05 & log2FoldChange > 0.58)#LFC > 1.5 or upReg in Pip
nrow(upReg)
head(upReg)

downReg <- subset(resLFC, padj < 0.05 & log2FoldChange < -0.58)#LFC < 1.5 or downReg in Pip
nrow(downReg)
head(downReg)

#Create plots based on the LFC, which minimizes noise from low read counts
pdf("plotMA_gravidVPipEvanF.pdf",width=6,height=6,paper='special')
plotMA(resLFC, ylim=c(-3,3))
dev.off()

pdf("EV_gravidVPipEvanF.pdf",width=8,height=6,paper='special')
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c('CPIJ003456'),
                selectLab = NA,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)
dev.off()

#Standard analysis - non LFC data 
res05 <- results(dds, alpha=0.05, contrast = c("condition", "PipEvanF", "CALgravidF"))
summary(res05)
sum(res05$padj < 0.05 & abs(resLFC$log2FoldChange) > 0.58, na.rm=TRUE)

#getting full gene list
full_genes_CalgravidF_PipEvanF <- data.frame(c(paste(resLFC@rownames, sep = "", "-RA")), resLFC@listData$log2FoldChange, resLFC@listData$padj)#paste adds VB format to IDs
write.table(full_genes_CalgravidF_PipEvanF, file = "CALgravidF_v_PipEvanF_AllGenes_LFCS.txt", 
            col.names = F, row.names = F, sep = "\t")

full_genes_CalgravidF_PipEvanF <- data.frame(c(paste(res05@rownames, sep = "", "-RA")), res05@listData$log2FoldChange, res05@listData$padj)#paste adds VB format to IDs
write.table(full_genes_CalgravidF_PipEvanF, file = "CALgravidF_v_PipEvanF_AllGenes_nonLFCS.txt", 
            col.names = F, row.names = F, sep = "\t")


#Save a list of significant DE genes
resSig <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) > 0.58)
write.table(data.frame(c(paste(resSig@rownames, sep = "", "-RA")), resSig@listData$log2FoldChange, resSig@listData$padj),
            file="CALgravidF_v_PipEvanF_LFCS_padj05.txt", row.names = F, col.names = F, sep = "\t") 

resSig <- subset(res05, padj < 0.05 & abs(log2FoldChange) > 0.58 )
write.table(data.frame(c(paste(resSig@rownames, sep = "", "-RA")), resSig@listData$log2FoldChange, resSig@listData$padj),
            file="CALGravidF_v_PipEvanF_nonLFCS_padj05.txt", row.names = F, col.names = F, sep = "\t") 


#VennDiag DEGs
dat_1 <- read.table("CALgravidF_v_CALparousF_LFCS_padj05.txt", header = F)
dat_2 <- read.table("CALparousF_v_PipEvanF_LFCS_padj05.txt", header = F)
dat_3 <- read.table("CALgravidF_v_PipEvanF_LFCS_padj05.txt", header = F)


l1 <- dat_1[, c(1,2)]
l2 <- dat_2[, c(1,2)]
l3 <- dat_3[, c(1,2)]

VennDiag <- GOVenn(l1,l2,l3, label=c('BG1 gravid v parous','AG2 v BG1 parous','AG2 v BG1 gravid'), plot = F)
print(VennDiag$plot)


######After looking at the first PCA, I wasn't sure about whether to drop BG_parous4.  In the end, I did not because it did not greatly impact the numbers of differentially expressed genes I recovered.
#Load in read counts and assign them sample labels
outputPrefix<-("Culex_DEseq_dropped")
sampleFiles<-c("M1-1_S1_htseq","M1-2_S2_htseq","M1-3_S3_htseq","M1-4_S4_htseq","M2-1_S5_htseq","M2-2_S6_htseq","M2-3_S7_htseq","M4-1_S13_htseq","M4-2_S14_htseq","M4-3_S15_htseq","M4-4_S16_htseq")
sampleNames<-c("BG_Gravid1", "BG_Gravid2", "BG_Gravid3", "BG_Gravid4", "BG_Parous1", "BG_Parous2", 
               "BG_Parous3", "AG1", "AG2", "AG3", "AG4")
sampleCondition<-c("CALgravidF",        "CALgravidF",   "CALgravidF",   "CALgravidF",   "CALparousF",   "CALparousF",   "CALparousF",   "PipEvanF",     "PipEvanF",     "PipEvanF",     "PipEvanF")
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)
treatments<-c("CALgravidF",     "CALparousF",   "PipEvanF")

#Create DESeq Data
ddsHTSeq_dropped <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                               directory = hiQual_Ex,
                                               design= ~ condition)

colData(ddsHTSeq_dropped)$condition <- factor(colData(ddsHTSeq_dropped)$condition, levels=treatments)
dim(ddsHTSeq_dropped)#getting num genes in dataset and verifying num samples


#Prefilter so only genes with atleast 10 reads in atleast 4 samples are considered.
keep <- rowSums(counts(ddsHTSeq_dropped) >= 10) >= 4
ddsHTSeq_dropped <- ddsHTSeq_dropped[keep,]
dim(ddsHTSeq_dropped)#new filtered num genes in dataset and verifying num samples

ddsHTSeq_dropped$condition <- factor(ddsHTSeq_dropped$condition, levels = treatments)


#new pca
rlogcounts_dropped <- rlog(counts(ddsHTSeq_dropped))

#run PCA
pcDat_dropped <- prcomp(t(rlogcounts_dropped))

#basic plot
autoplot(pcDat_dropped)

#plot for pub
pdf("PCA_Treatments_dropped.pdf",width=6,height=6,paper='special')
autoplot(pcDat_dropped,
         data = ddsHTSeq_dropped$colData, 
         colour=as.numeric(factor(ddsHTSeq_dropped$condition)), 
         shape=FALSE,
         label.size=6, xlim = c(-0.4, 0.5)) + theme_bw()
dev.off()


####This is the DeSeq2 analysis BG_gravid v BG_parous without the BG_parous4
dds_dropped<-DESeq(ddsHTSeq_dropped)
res_dropped<-results(dds_dropped)

###Gravid F vs. Parous F--Use LFC for gene ranking and visualization, and use the p-values from the non-LFC
resLFC_dropped <- lfcShrink(dds_dropped, contrast = c("condition", "CALgravidF", "CALparousF"), type="ashr")

#Look at summary values
summary(resLFC_dropped)

#how many significantly DE genes? The default p-value cutoff is 0.1
sum(resLFC_dropped$padj < 0.1 & abs(resLFC_dropped$log2FoldChange) > 0.58, na.rm=TRUE)
sum(resLFC_dropped$padj < 0.05, na.rm=TRUE)
sum(resLFC_dropped$padj < 0.05 & abs(resLFC_dropped$log2FoldChange) > 0.58, na.rm=TRUE)

#Look at genes when changing the p-value cutoff--Use p-values from the non LFC data 
res05_dropped <- results(dds_dropped, alpha=0.05, contrast = c("condition", "CALgravidF", "CALparousF"))
summary(res05_dropped)

sum(res05_dropped$padj < 0.05 & abs(resLFC_dropped$log2FoldChange) > 0.58, na.rm=TRUE)


#Create plots based on the LFC, which minimizes noise from low read counts
pdf("plotMA_gravidVparous_dropped.pdf",width=6,height=6,paper='special')
plotMA(resLFC_dropped, ylim=c(-3,3))
dev.off()

#Save a list of significant DE genes
resSig_dropped <- subset(res05_dropped, padj < 0.05  & abs(log2FoldChange) > 0.58)
write.csv(as.data.frame(resSig_dropped), 
          file="CALgravidF_v_CALparousF_dropped_nonLFCS_padj05.csv")

resSig_dropped <- subset(resLFC_dropped, padj < 0.05 & abs(log2FoldChange) > 0.58)
write.csv(as.data.frame(resSig_dropped), 
          file="CALgravidF_v_CALparousF_dropped_LFCS_padj05.csv")

pdf("EV_gravidVparous_dropped.pdf",width=8,height=6,paper='special')
EnhancedVolcano(resLFC_dropped,
                lab = rownames(resLFC_dropped),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c('CPIJ003456'),
                selectLab = NA,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)
dev.off()

######Here rerunning DeSeq2 for Parous v pip without BG_parous 4 outlier

resLFC_dropped <- lfcShrink(dds_dropped, contrast = c("condition", "PipEvanF", "CALparousF"), type="ashr")

#Look at summary values
summary(resLFC_dropped)

#how many significantly DE genes? The default p-value cutoff is 0.1
sum(resLFC_dropped$padj < 0.1  & abs(resLFC_dropped$log2FoldChange) > 0.58, na.rm=TRUE)
sum(resLFC_dropped$padj < 0.05, na.rm=TRUE)
sum(resLFC_dropped$padj < 0.05 & abs(resLFC_dropped$log2FoldChange) > 0.58, na.rm=TRUE)

#Look at genes when changing the p-value cutoff--Use p-values from the non LFC data 
res05_dropped <- results(dds_dropped, alpha=0.05, contrast = c("condition", "PipEvanF", "CALparousF"))
summary(res05_dropped)
sum(res05_dropped$padj < 0.05  & abs(resLFC_dropped$log2FoldChange) > 0.58, na.rm=TRUE)

#Create plots based on the LFC, which minimizes noise from low read counts
pdf("plotMA_parousVpip_dropped.pdf",width=6,height=6,paper='special')
plotMA(resLFC_dropped, ylim=c(-3,3))
dev.off()

pdf("EV_parousVpip_dropped.pdf",width=8,height=6,paper='special')
EnhancedVolcano(resLFC_dropped,
                lab = rownames(resLFC_dropped),
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = c('CPIJ003456'),
                selectLab = NA,
                xlim = c(-10, 10),
                ylim = c(0,200),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)
dev.off()

#Save a list of significant DE genes
resSig_dropped <- subset(res05_dropped, padj < 0.05 & abs(log2FoldChange) > 0.58)
write.csv(as.data.frame(resSig_dropped), 
          file="CALparousF_v_PipEvanF_dropped_nonLFCS_padj05.csv")

resSig_dropped <- subset(resLFC_dropped, padj < 0.05  & abs(log2FoldChange) > 0.58)
write.csv(as.data.frame(resSig_dropped), 
          file="CALparousF_v_PipEvanF_dropped_LFC_padj05.csv")


