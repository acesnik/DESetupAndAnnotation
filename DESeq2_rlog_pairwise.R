options(repos = c(CRAN = "https://cran.rstudio.com"))
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
#install.packages("optparse")
library(optparse)
library("RColorBrewer")
library("ggplot2")

#Take in filenames
option_list = list(
  make_option(c("-i", "--input.dataframe"), type="character", default = NULL, 
              help="Dataset of estimated counts from RSEM in this case. Each column should correspond to the descending rows in the sample list."),
  make_option(c("-s", "--sample.list"), type="character", default = NULL, 
              help="Sample list. Each column can be a variable in the design. All entries must start with a letter."),
  make_option(c("-o", "--out.folder"), type="character", default=".", 
              help="output file name [default= %default]")
  make_option(c("-e", "--main.effect"), type="character", default="main_effect",
              help="main effect, prefix for filenames")
); 
opt = parse_args(OptionParser(option_list=option_list));

sample.data.file <- opt$input.dataframe
sample.list <- opt$sample.list
output.folder <- opt$out.folder
main.effect <- opt$main.effect

#Count matrix input
sample.data<-round(read.delim(sample.data.file, row.names=1, header = TRUE))
colnames(sample.data) <- NULL
sample.list<-read.delim(sample.list, row.names=1)
dds <- DESeqDataSetFromMatrix(countData = sample.data, colData = sample.list, design=~class)
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
distsRL <- dist(t(assay(rlt)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(de), paste(patient, type, class, sep=" : "))
pdf(paste(output.folder, "/", main.effect, "_", this.plot, ".pdf", sep=""), width = 6, height = 6)
heatmap(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()
this.plot <- "dispersionDiagnostic"
pdf(paste(output.folder,"/", main.effect, "_", this.plot, ".pdf", sep=""), width = 11, height = 6)
plotDispEsts(de)
dev.off()

this.plot <- "effectOnNormal"
res<-results(de, contrast=c("class","aggressive","indolent"))
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
