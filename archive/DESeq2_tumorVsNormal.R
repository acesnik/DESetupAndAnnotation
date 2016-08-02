source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

#Take in filenames
sample.data.file <- "/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/_paired_tumorVsNormal/rsemStar_mRNA_canon/isoform_mRNA_dataframe.txt"
sample.list <- "prostate_sampleList.txt"
output.folder <- "/Volumes/KemiTegel/ProstateTumorAnalysis/DataAnalysisR/_basic_aggVsInd_tumorNormalInteraction/rsemStar_mRNA_canon_vsd"
experiment.name <- "aggVsInd"

#Count matrix input
sample.data<-round(read.delim(sample.data.file, row.names=1, header = TRUE))
sample.list<-read.delim(sample.list, row.names=1)
dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design=~type+class+type:class)
de <- DESeq(dds)

#Diagnostic MA plot
pdf(paste(output.folder,"/",experiment.name,"tumorVsNormal_MAPlot.pdf", sep=""), width = 11, height = 6)
plotMA(de)
dev.off()

#Transcript-level hatmaps
rlt <- rlogTransformation(de, blind=TRUE)
vsd <- varianceStabilizingTransformation(de, blind=TRUE)

library("RColorBrewer")
library("ggplot2")
select <- order(rowMeans(counts(de, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf(paste(output.folder, "/", experiment.name, "_transcriptLevelHeatmap.pdf", sep=""), width = 11, height = 6)
heatmap(counts(de,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.off()

pdf(paste(output.folder, "/", experiment.name, "_transcriptLevelHeatmap_rLogTransformNormalization.pdf", sep=""), width = 11, height = 6)
heatmap(assay(rlt)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

pdf(paste(output.folder, "/", experiment.name, "_transcriptLevelHeatmap_varStabilizingTransform.pdf", sep=""), width = 11, height = 6)
heatmap(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

#Samplewise heatmap
distsRL <- dist(t(counts(de,normalized=TRUE)[select,]))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(patient, type, class, sep=" : "))
pdf(paste(output.folder, "/", experiment.name, "_samplewiseHeatmap.pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(patient, type, class, sep=" : "))
pdf(paste(output.folder, "/", experiment.name, "_samplewiseHeatmap.pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(patient, type, class, sep=" : "))
pdf(paste(output.folder, "/", experiment.name, "_samplewiseHeatmap.pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#Dispersion diagnostic plot
pdf(paste(output.folder,"/tumorVsNormal_dispersionDiagnostic.pdf", sep=""), width = 11, height = 6)
plotDispEsts(de)
dev.off()

#Export results
res<-results(de, contrast=c("class","aggressive","indolent"))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df), file=paste(output.folder,"/aggVsInd_normal.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res), file=paste(output.folder,"/aggVsInd_normal_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

res<-results(de, list(c("class_indolent_vs_aggressive","typetumor.classindolent")))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df), file=paste(output.folder,"/aggVsInd_tumor.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res), file=paste(output.folder,"/aggVsInd_tumor_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

res<-results(de, name="typetumor.classindolent")
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df), file=paste(output.folder,"/aggVsInd_differentEffect.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
