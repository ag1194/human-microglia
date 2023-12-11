##Code to do DESeq2 analysis from read counts
##Written by Aysegul Guvenek
##November 2023

library("DESeq2")
library("GenomicAlignments")
library("GenomicRanges")
library("rtracklayer")
library("Rsamtools")
library("GenomicFeatures")
library("dplyr")
library("ggplot2")

rm(list = ls())
par(ask=TRUE)

lc=12 #number of control samples
lt=12 #number of treatment samples

data = read.table("all.Gene_reads_TPM.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
data_subset = data[,c("gene_symbol",grep("CDSreads", colnames(data),value = TRUE))]

data_subset00 = data_subset[,c(grep("cntrl", colnames(data_subset),value = TRUE))]
data_subset01 = data_subset[, !grepl("cntrl", names(data_subset))]
rownames(data_subset01)= data_subset01$gene_symbol
rownames(data_subset00)= data_subset01$gene_symbol
data_subset01<- data_subset01[-1]
sample_c<- colnames(data_subset00)
sample_t<- colnames(data_subset01)
sample_c1 <-sub("_CDSreads","", sample_c)
sample_t1<-sub("_CDSreads","", sample_t)

dfinput0<-cbind(data_subset00, data_subset01)
dfinput0 $min = as.numeric (apply(dfinput0, 1, FUN= min))
dfinput0 = subset(dfinput0, min>0)
dfinput0 = dfinput0[-length(dfinput0)]


##DESeq2 analysis

condition = factor(c((rep("control",lc)),(rep("treated",lt))))
sampleTable = data.frame(row.names = c(colnames(data_subset00), colnames(data_subset01)),
                         condition = condition)
dds <- DESeqDataSetFromMatrix(countData = dfinput0,
                              colData = sampleTable,
                              design= ~  condition)
dds <- DESeq(dds)
res <- results(dds)
result_dds<- as.data.frame(res)
result_dds2<-result_dds[,c("pvalue","padj","log2FoldChange")]

dfinput2=as.data.frame(counts(dds, normalized=TRUE))
dfinput2<- cbind(dfinput2, result_dds2)
dfinput2 $change.padj ="NC"
dfinput2 $change.padj=ifelse((dfinput2 $padj <0.05 & (dfinput2 $log2FoldChange>1)), "UP", dfinput2 $change.padj)
dfinput2 $change.padj=ifelse((dfinput2 $padj <0.05 & (dfinput2 $log2FoldChange< -1)), "DN", dfinput2 $change.padj)

#####plotting

res2 = result_dds
res3<-subset(res2, log2FoldChange !="NaN" & log2FoldChange !="Inf" & log2FoldChange !="-Inf")

res3$change.padj ="NC"
res3$change.padj=ifelse((res3$padj <0.05 & (res3$log2FoldChange>1)), "UP", res3$change.padj)
res3$change.padj=ifelse((res3$padj <0.05 & (res3$log2FoldChange< -1)), "DN", res3$change.padj)

res3$color="black"
res3$color[res3$change.padj=="UP"]="red"
res3$color[res3$change.padj=="DN"]="blue"

blue=table(res3 $change.padj)[[1]]
gray=table(res3 $change.padj)[[2]]
red=table(res3 $change.padj)[[3]]

title=table(res3$color)
table(res3$color)

xi="Microglia_rev"
zi="Parietal-cortex"

pdf(paste0(xi,'_vs_',zi,".DESeq.Galatro.pdf"),width=8, height=8)
with(res3, plot(log2FoldChange, -log10(padj), main=paste("DeSeq FC:2 Padj:0.05 ", xi," vs ", zi, "\n red: ", red," ", " gray: ", gray, " blue: ", blue), col=color, pch = ifelse((color =="black"), 20,19),  xlab="log2 Fold Change", ylab="-log10 padj"))
dev.off()

write.table(dfinput2, paste0(xi,'_vs_',zi,".all.CDS_DESeq.txt"), sep="\t",row.names=T,quote=FALSE)


