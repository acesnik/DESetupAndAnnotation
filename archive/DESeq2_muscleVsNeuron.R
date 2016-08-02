#source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

#Count matrix input
worm.data<-round(read.delim("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/isoform_dataframe_noCtrl.txt", row.names=1))
sample.list<-read.delim("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/160525_scRNAseq_SAMPLE_LIST_noCtrl.txt", row.names=1)
dds <- DESeqDataSetFromMatrix(countData = worm.data, colData = sample.list, design = ~cell_type)
de <- DESeq(dds)

#Diagnostic MA plot
pdf("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron_MAPlot.pdf", width = 11, height = 6)
plotMA(de)
dev.off()

#Transcript-level hatmaps
rlt <- rlogTransformation(de, blind=TRUE)
vsd <- varianceStabilizingTransformation(de, blind=TRUE)

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(de, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron_transcriptLevelHeatmap.pdf", width = 11, height = 6)
heatmap(counts(de,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.off()
pdf("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron_transcriptLevelHeatmap_rLogTransformNormalization.pdf", width = 11, height = 6)
heatmap(assay(rlt)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
pdf("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron_transcriptLevelHeatmap_varStabilizingTransform.pdf", width = 11, height = 6)
heatmap(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.off()

#Samplewise heatmap
distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de),
                                       paste(cell_type, cell_function, sep=" : "))
pdf("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron_samplewiseHeatmap.pdf", width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#Dispersion diagnostic plot
pdf("/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron_dispersionDiagnostic.pdf", width = 11, height = 6)
plotDispEsts(de)
dev.off()

#Export results
res<-results(de)
write.table(res, file="/Volumes/KemiTegel/WormSingleCell/160615_RsemStar_OnlyWormRef/MuscleVsNeuron_deseq2/muscleVsNeuron.txt", 
            sep="\t", quote=FALSE)
