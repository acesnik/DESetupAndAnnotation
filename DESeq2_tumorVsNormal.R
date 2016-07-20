source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

#Take in filenames
args = commandArgs(trailingOnly=TRUE)
usage = "Rscript DESeq_tumorVsNormal.R comparison dataframe.txt sample_list.txt output_folder\n
comparison can be ~patient; ~tissue_type, ~tissue_class, or any combination of those."
if (length(args) == 0) { stop(usage, call.=FALSE)
comparison <- args[1]
sample.data.file <- args[2]
sample.list <- args[3]
output.folder <- args[4] 

#Count matrix input
sample.data<-round(read.delim(sample.data.file, row.names=1, header = TRUE))
sample.list<-read.delim(sample.list, row.names=1)
if (comparison == "~tissue_type") dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~tissue_type)
else if (comparison == "~tissue_class") dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~tissue_class)
else if (comparison == "~patient")  dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~patient)
else if (comparison == "~patient+tissue_type") dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~patient+tissue_type)
else if (comparison == "~patient+tissue_class")  dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~patient+tissue_class)
else if (comparison == "~patient+tissue_class+tissue+type")  dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design = ~patient+tissue_class+tissue+type)
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
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df), file=paste(output.folder,"/tumorVsNormal.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res), file=paste(output.folder,"/tumorVsNormal_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)