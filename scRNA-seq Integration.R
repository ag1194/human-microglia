#####INTEGRATE Seurat objects

##################################################
##########to increase the capacity###############
library(future)
# check the current active plan
plan()
# change the current plan to access parallelization
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 12000 * 2048^32)
plan()
##################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)



####### upload individually analyzed scRNA-seq datasets###### Each are SCTransform by Seurat

msNew <- readRDS("/mnt/efs_v2/tag_neuro/users/aysegul.guvenek/analyze/Microglia/human/Schirmer_Rowitch/NEwmsSCT.rds", refhook = NULL)
ad <- readRDS("/mnt/efs_v2/tag_neuro/users/aysegul.guvenek/analyze/Microglia/human/integrated/Feb2022/adSCT.rds", refhook = NULL)
ms <- readRDS("/mnt/efs_v2/tag_neuro/users/aysegul.guvenek/analyze/Microglia/human/Parkinson_GSE157783/msSCT.rds", refhook = NULL)
prk <- readRDS("/mnt/efs_v2/tag_neuro/users/aysegul.guvenek/analyze/Microglia/human/Parkinson_GSE157783/parkinsonSCT.rds", refhook = NULL)


######### Integrate all 4
reference.list <-  list(ad, ms, prk, msNew)
msadprk.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = msadprk.features, 
    verbose = FALSE)
msadprk.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", 
    anchor.features = msadprk.features, verbose = FALSE)
msadprk.integrated <- IntegrateData(anchorset = msadprk.anchors, normalization.method = "SCT", verbose = FALSE)
 
saveRDS(msadprk.integrated, file = "integrated_4datasets.rds")



all.integrated <- msadprk.integrated
all.integrated <-subset(all.integrated, CSF1R > 1 | BIN1>1 | P2RY12 >1) #### to subset microglia
all.integrated <- RunPCA(all.integrated, verbose = FALSE)
all.integrated <- FindNeighbors(all.integrated, dims = 1:10)
#all.integrated <- FindClusters(all.integrated, resolution = 0.6)
all.integrated <- FindClusters(all.integrated, resolution = 0.4)

all.integrated <- RunUMAP(all.integrated, dims = 1:30)
#all.integrated <- RunUMAP(all.integrated, dims = 1:10)


#######plotting#######
#### UMAP plots  #####

p2 <- DimPlot(all.integrated, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE)

pdf(paste0("./IntegratedSCt_3k_MSADNewMSPrk.pdf"),width=8, height=8)
p2 & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    override.aes = list(size = 3)))
dev.off()


p2 <- DimPlot(all.integrated, reduction = "umap", group.by = "ident",  repel = TRUE)

pdf(paste0("./IntegratedSCt_3k_MSADNewMSPrk1.pdf"),width=8, height=8)
p2 & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    override.aes = list(size = 3)))
dev.off()

plots <- DimPlot(all.integrated, group.by = "orig_id2")


plots <- DimPlot(all.integrated, group.by = "orig_id2", cols = c( "brown", "darkgoldenrod1", "darkorange","darkgoldenrod1","darkgoldenrod1","coral","darkgoldenrod1","green","aquamarine","deepskyblue2","dodgerblue","black","darkseagreen1","red","red","red","red","magenta"))


pdf(paste0("./IntegratedSCt_3k_MSADNewMSPrk_origID.pdf"),width=8, height=8)
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    override.aes = list(size = 3)))
dev.off()



######### Find Markers and make heatmap ################

pbmc.markers <- FindAllMarkers(all.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

write.csv(pbmc.markers, "MSADPrkMsNew_allClusters.markers")

library(ggplot2)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf(paste0("MSADPrkMsNew-clusterGenes_heatmap.markers.pdf"),width=12, height=16)
DoHeatmap(all.integrated, features = top10$gene)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, indent.2=c(0,1,2,3), min.pct = 0.25)

write.csv(cluster4, "CDAM.markers.csv")

########################################
