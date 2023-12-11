
library(pheatmap)

heatmapgenes<- read.csv(ordered_c,"upset_FCby1.5_ordered.csv")
heatmapgenes1 <- na.omit(heatmapgenes)

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max1/paletteLength, 4, length.out=floor(paletteLength/2)))


heatmapgenes2 <- heatmapgenes1[2:50]
rownames(heatmapgenes2)= heatmapgenes $gene

l1<- length(heatmapgenes2[[1]])
pdf('combined_heatmap_woNA_upsetbased_proteincodingOnly.pdf',width=40, height=20)
pheatmap(t(heatmapgenes2),cluster_rows =F, cluster_cols =F, scale="none", fontsize = 15, show_rownames=T, show_colnames=F, dpi=450, width = 50, height = 20, ,color=myColor, breaks=myBreaks, treeheight_row=50, treeheight_col=300,main=paste("log2FC normalized within Dataset,mRNA genes p<0.05, fc 1.5, noNA", l1))
dev.off()


