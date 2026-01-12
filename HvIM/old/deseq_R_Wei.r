#T4
readscount <- read_excel('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Input/raw data_readcount_no rep2.xlsx', sheet = "T4")
colData <- read_excel('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Input/colData_no rep2.xlsx', sheet = "T4")
row.names(readscount) <- geneIDs$geneID #Assigning row names from geneIDs
condition <-factor(c(c("Control", "Flag22", "Pnic", "Flag22+Pnic"), c("Control", "Flag22", "Pnic"), c("Control", "Flag22", "Pnic", "Flag22+Pnic")))
replicate <- factor(c(rep("Rep1", 4), rep("Rep3", 3), rep("Rep0", 4)))
colData
head(readscount)
condition
replicate
dds <- DESeqDataSetFromMatrix (readscount, colData, design = ~replicate + condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

#PCA analysis###
pdf("T4_PCA_no rep2.pdf", width = 6, height = 6)
vsdata <- vst(dds, blind=FALSE) #variance stabilizing transformation
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata),vsdata$condition)
plotPCA(vsdata, intgroup = "replicate")

vsdata <- vst(dds, blind=FALSE)
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata),vsdata$replicate)
plotPCA(vsdata, intgroup = "condition")
dev.off()


#Plot of expression values
readscount_new = assay(vsdata)
head(readscount_new)
par(cex = 0.7)
n.sample=ncol(readscount)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(readscount, col = cols, main="expression value",las=2, labels = FALSE)
boxplot(readscount_new, col = cols,main="expression value",las=2, labels = FALSE)
hist(readscount_new, main = "Readscount Histgram")

# Differential expression analysis between each treatment and control
dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
dds_norm$condition   #make sure levels are treat/untreat，otherwise needs to explicitely say it when define results
res_Pnic <- results(dds_norm, contrast = c("condition","Pnic","Control"), cooksCutoff = FALSE) #alpha=0.05 for padj; cookCutoff (for filter out outlier)= false cuz there are lots of outlier
dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
dds_norm$condition   #make sure levels are treat/untreat，otherwise needs to explicitely say it when define results
res_Flag22 <- results(dds_norm, contrast = c("condition","Flag22","Control"), cooksCutoff = FALSE) #alpha=0.05 for padj; cookCutoff (for filter out outlier)= false cuz there are lots of outlier
dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
dds_norm$condition   #make sure levels are treat/untreat，otherwise needs to explicitely say it when define results
res_Flag22_Pnic <- results(dds_norm, contrast = c("condition","Flag22+Pnic","Control"), cooksCutoff = FALSE) #alpha=0.05 for padj; cookCutoff (for filter out outlier)= false cuz there are lots of outlier

summary(res_Pnic)
summary(res_Flag22)
summary(res_Flag22_Pnic)

#Filter out the ones with significance with padj < 0.05
sum(res_Pnic$padj<0.05, na.rm = TRUE)
sum(res_Flag22$padj<0.05, na.rm = TRUE)
sum(res_Flag22_Pnic$padj<0.05, na.rm = TRUE)

res_Pnic_data <- merge(as.data.frame(res_Pnic),
                       as.data.frame(counts(dds_norm,normalize=TRUE)),
                       by="row.names",sort=FALSE)
res_Flag22_data <- merge(as.data.frame(res_Flag22),
                         as.data.frame(counts(dds_norm,normalize=TRUE)),
                         by="row.names",sort=FALSE)
res_Flag22_Pnic_data <- merge(as.data.frame(res_Flag22_Pnic),
                              as.data.frame(counts(dds_norm,normalize=TRUE)),
                              by="row.names",sort=FALSE)

up_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22_PnicDEG <- subset(res_Flag22_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22_PnicDEG <- subset(res_Flag22_Pnic_data, padj < 0.05 & log2FoldChange < -1)


write.csv(res_Pnic_data, "T4_all_Pnic_no rep2.csv") #Not filtered data, for volcano
write.csv(up_PnicDEG, "T4_up_Pnic_no rep2.csv")
write.csv(down_PnicDEG, "T4_down_Pnic_no rep2.csv")

write.csv(res_Flag22_data, "T4_all_Flag22_no rep2.csv") #Not filtered data, for volcano
write.csv(up_Flag22DEG, "T4_up_Flag22_no rep2.csv")
write.csv(down_Flag22DEG, "T4_down_Flag22_no rep2.csv")

write.csv(res_Flag22_Pnic_data, "T4_all_Flag22_Pnic_no rep2.csv") #Not filtered data, for volcano
write.csv(up_Flag22_PnicDEG, "T4_up_Flag22_Pnic_no rep2.csv")
write.csv(down_Flag22_PnicDEG, "T4_down_Flag22_Pnic_no rep2.csv")

