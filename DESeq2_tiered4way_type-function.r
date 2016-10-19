#Take in filenames
options(repos = c(CRAN = "https://cran.rstudio.com"))
#install.packages("optparse")
library(optparse)
library("RColorBrewer")
library("ggplot2")
option_list = list(
  make_option(c("-i", "--input.dataframe"), type="character", default = NULL, 
              help="Dataset of estimated counts from RSEM in this case. Each column should correspond to the descending rows in the sample list."),
  make_option(c("-s", "--sample.list"), type="character", default = NULL, 
              help="Sample list. Each column can be a variable in the design. All entries must start with a letter."),
  make_option(c("-o", "--out.folder"), type="character", default=".", 
              help="output file name [default= %default]")
); 
opt = parse_args(OptionParser(option_list=option_list));

sample.data.file <- opt$input.dataframe
sample.list <- opt$sample.list
output.folder <- opt$out.folder

### MAIN EFFECT: CELL TYPE AND FUNCTION ###
# note: remember to change the design, too
# note: remember to change the heatmap labels, too
main.effect <- "cellTypeAndFunction"
comparison1 <- "cell"
up1 <- "NeuronExcitatory"
down1 <- "NeuronInhibitory"
up2 <- "NeuronExcitatory"
down2 <- "Muscle"
up3 <- "NeuronInhibitory"
down3 <- "Muscle"
up4 <- "NeuronAVK"
down4 <- "Muscle"

#Count matrix input
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
sample.data<-round(read.delim(sample.data.file, row.names=1, header = TRUE))
colnames(sample.data) <- NULL
sample.list<-read.delim(sample.list, row.names=1)
dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design=~cell)
de <- DESeq(dds)

this.plot <- "diagnosticMAplot"
pdf(paste(output.folder,"/",main.effect,"_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotMA(de)
dev.off()

# Regularized logarithmic transformation, which does sample-wise normalization of counts
rlt <- rlogTransformation(de, blind=TRUE)
save.image(file=paste(output.folder, "/", main.effect, ".RData", sep=""))

this.plot <- "rlt_transcriptLevelHeatmap"
select <- order(rowMeans(counts(de, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
heatmap(assay(rlt)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
this.plot <- "rlt_samplewiseHeatmap"
de@colData$sample <- rownames(sample.list)
distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(sample, cell, sep=" : "))
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()
this.plot <- "dispersionDiagnostic"
pdf(paste(output.folder,"/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotDispEsts(de)
dev.off()

this.plot <- paste(up1, "_vs_", down1)
res<-results(de, contrast=c(comparison1, up1, down1))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
res_df$negLogPadj <- -log(res_df$padj)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, "_volcano.pdf", sep=""))
ggplot(res_df, aes(y=negLogPadj, x=log2FoldChange)) + labs(x="Log-2 Fold Change", y="-log(p-value, adjusted)") + geom_point()
dev.off()
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df),
file=paste(output.folder, "/", main.effect, "_", this.plot, ".txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res),
file=paste(output.folder, "/", main.effect, "_", this.plot, "_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

this.plot <- paste(up2, "_vs_", down2)
res<-results(de, contrast=c(comparison1, up2 , down2))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
res_df$negLogPadj <- -log(res_df$padj)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, "_volcano.pdf", sep=""))
ggplot(res_df, aes(y=negLogPadj, x=log2FoldChange)) + labs(x="Log-2 Fold Change", y="-log(p-value, adjusted)") + geom_point()
dev.off()
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df),
file=paste(output.folder, "/", main.effect, "_", this.plot, ".txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res),
file=paste(output.folder, "/", main.effect, "_", this.plot, "_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

this.plot <- paste(up3, "_vs_", down3)
res<-results(de, contrast=c(comparison1, up3 , down3))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
res_df$negLogPadj <- -log(res_df$padj)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, "_volcano.pdf", sep=""))
ggplot(res_df, aes(y=negLogPadj, x=log2FoldChange)) + labs(x="Log-2 Fold Change", y="-log(p-value, adjusted)") + geom_point()
dev.off()
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df),
file=paste(output.folder, "/", main.effect, "_", this.plot, ".txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res),
file=paste(output.folder, "/", main.effect, "_", this.plot, "_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

this.plot <- paste(up4, "_vs_", down4)
res<-results(de, contrast=c(comparison1, up4 , down4))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
res_df$negLogPadj <- -log(res_df$padj)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, "_volcano.pdf", sep=""))
ggplot(res_df, aes(y=negLogPadj, x=log2FoldChange)) + labs(x="Log-2 Fold Change", y="-log(p-value, adjusted)") + geom_point()
dev.off()
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df),
            file=paste(output.folder, "/", main.effect, "_", this.plot, ".txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res),
            file=paste(output.folder, "/", main.effect, "_", this.plot, "_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

save.image(file=paste(output.folder, "/", main.effect, ".RData", sep=""))

### MAIN EFFECT: CELL TYPE ###
main.effect <- "cellType"
comparison <- "cell.type"
up <- "Neuron"
down <- "Muscle"

dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design=~cell.type)
de <- DESeq(dds)
this.plot <- "diagnosticMAplot"
pdf(paste(output.folder,"/",main.effect,"_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotMA(de)
dev.off()

# Regularized logarithmic transformation, which does sample-wise normalization of counts
rlt <- rlogTransformation(de, blind=TRUE)
save.image(file=paste(output.folder, "/", main.effect, ".RData", sep=""))


this.plot <- "rlt_transcriptLevelHeatmap"
select <- order(rowMeans(counts(de, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
heatmap(assay(rlt)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
this.plot <- "rlt_samplewiseHeatmap"
de@colData$sample <- rownames(sample.list)
distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(sample, cell, sep=" : "))
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()
this.plot <- "dispersionDiagnostic"
pdf(paste(output.folder,"/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotDispEsts(de)
dev.off()

this.plot <- paste(up, "_vs_", down)
res<-results(de, contrast=c(comparison, up, down))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
res_df$negLogPadj <- -log(res_df$padj)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, "_volcano.pdf", sep=""))
ggplot(res_df, aes(y=negLogPadj, x=log2FoldChange)) + labs(x="Log-2 Fold Change", y="-log(p-value, adjusted)") + geom_point()
dev.off()
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df),
file=paste(output.folder, "/", main.effect, "_", this.plot, ".txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res),
file=paste(output.folder, "/", main.effect, "_", this.plot, "_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

save.image(file=paste(output.folder, "/", main.effect, ".RData", sep=""))

### MAIN EFFECT: CELL FUNCTION AND PROJECTION ###
main.effect <- "cellFunctionAndProjection"
comparison2 <- "function.projection"
up5 <- "ExcitatoryAnterior"
down5 <- "ExcitatoryPosterior"

dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design=~function.projection)
de <- DESeq(dds)
this.plot <- "diagnosticMAplot"
pdf(paste(output.folder,"/",main.effect,"_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotMA(de)
dev.off()

# Regularized logarithmic transformation, which does sample-wise normalization of counts
rlt <- rlogTransformation(de, blind=TRUE)
save.image(file=paste(output.folder, "/", main.effect, ".RData", sep=""))


this.plot <- "rlt_transcriptLevelHeatmap"
select <- order(rowMeans(counts(de, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
heatmap(assay(rlt)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
dev.off()
this.plot <- "rlt_samplewiseHeatmap"
de@colData$sample <- rownames(sample.list)
distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(sample, cell, function.projection, sep=" : "))
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()
this.plot <- "dispersionDiagnostic"
pdf(paste(output.folder,"/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotDispEsts(de)
dev.off()

this.plot <- paste(up5, "_vs_", down5)
res<-results(de, contrast=c(comparison2, up5, down5))
res_df <- data.frame(res@listData, row.names = res@rownames)
res_df <- res_df[order(res_df$padj),]
res_df$negLogPadj <- -log(res_df$padj)
pdf(paste(output.folder, "/", main.effect, "_", this.plot, "_volcano.pdf", sep=""))
ggplot(res_df, aes(y=negLogPadj, x=log2FoldChange)) + labs(x="Log-2 Fold Change", y="-log(p-value, adjusted)") + geom_point()
dev.off()
top_res <- subset.data.frame(res_df, res_df$padj <= 0.01)
top_res <- top_res[order(-abs(top_res$log2FoldChange)),]
write.table(data.frame("transcript_id"=rownames(res_df), res_df),
            file=paste(output.folder, "/", main.effect, "_", this.plot, ".txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(data.frame("transcript_id"=rownames(top_res), top_res),
            file=paste(output.folder, "/", main.effect, "_", this.plot, "_topResults.txt", sep=""), row.names=FALSE, sep="\t", quote=FALSE)

save.image(file=paste(output.folder, "/", main.effect, ".RData", sep=""))

q()
