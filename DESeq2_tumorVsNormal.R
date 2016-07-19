#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library(DESeq2)

#Take in filenames
sample.list <- '/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/rsemStar_mRNA_cuff/prostate_sub_sampleList.txt'
sample.data.file <- '/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/rsemStar_rRNAd_cuff/isoform_dataframe.txt'
output.folder <- '/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/rsemStar_rRNAd_cuff'

#Count matrix input
sample.data<-round(read.delim(sample.data.file, row.names=1, header = TRUE))
sample.list<-read.delim(sample.list, row.names=1)
dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~tissue_type)
de <- DESeq(dds)

#Diagnostic MA plot
pdf(paste(output.folder,"/tumorVsNormal_MAPlot.pdf", sep=""), width = 11, height = 6)
plotMA(de)
dev.off()

#Transcript-level hatmaps
rlt <- rlogTransformation(de, blind=TRUE)
vsd <- varianceStabilizingTransformation(de, blind=TRUE)

library("RColorBrewer")
library("ggplot2")
select <- order(rowMeans(counts(de, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf(paste(output.folder,"/tumorVsNormal_transcriptLevelHeatmap.pdf", sep=""), width = 11, height = 6)
heatmap(counts(de,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.off()
pdf(paste(output.folder,"/tumorVsNormal_transcriptLevelHeatmap_rLogTransformNormalization.pdf", sep=""), width = 11, height = 6)
heatmap(assay(rlt)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
pdf(paste(output.folder,"/tumorVsNormal_transcriptLevelHeatmap_varStabilizingTransform.pdf", sep=""), width = 11, height = 6)
heatmap(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

#Samplewise heatmap
distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(patient, tissue_type, tissue_class, sep=" : "))
pdf(paste(output.folder,"/tumorVsNormal_samplewiseHeatmap.pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#Dispersion diagnostic plot
pdf(paste(output.folder,"/tumorVsNormal_dispersionDiagnostic.pdf", sep=""), width = 11, height = 6)
plotDispEsts(de)
dev.off()

#Export results
res<-results(de)
write.table(res, file=paste(output.folder,"/tumorVsNormal.txt", sep=""), 
            sep="\t", quote=FALSE)
