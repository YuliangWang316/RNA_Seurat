library(patchwork)
library(Seurat)
library(dplyr)
library(future)

pbmc<-readRDS("pbmc.rds")
pbmc_lda<-readRDS("d:/20231201/pbmc_lda.rds")

Idents(pbmc_lda)<-pbmc_lda$Myc
pbmc_Myc.markers<-FindAllMarkers(pbmc_lda,only.pos = TRUE,min.pct = 0)
pbmc_Myc.markers$diff <- pbmc_Myc.markers$pct.1 - pbmc_Myc.markers$pct.2
pbmc_Myc.markers_new<-pbmc_Myc.markers[which(pbmc_Myc.markers$cluster == "Myc" & pbmc_Myc.markers$p_val_adj < 0.05 & pbmc_Myc.markers$pct.1> 0.7 & pbmc_Myc.markers$diff > 0.3),]
Idents(pbmc_lda)<-pbmc_lda$Irf4
pbmc_Irf4.markers<-FindAllMarkers(pbmc_lda,only.pos = TRUE,min.pct = 0)
Idents(pbmc_lda)<-pbmc_lda$Ly75
pbmc_Ly75.markers<-FindAllMarkers(pbmc_lda,only.pos = TRUE,min.pct = 0)
Idents(pbmc_lda)<-pbmc_lda$Kdm6b
pbmc_Kdm6b.markers<-FindAllMarkers(pbmc_lda,only.pos = TRUE,min.pct = 0)

pbmc_ICA<-RunICA(pbmc,features = )
pbmc_PCA<-RunPCA(pbmc,features = )
pbmc_Neighbor_PCA <-FindNeighbors(pbmc_PCA,features = ,reduction = "pca",k.param = 10,dims = 1:4)
pbmc_Neighbor_PCA_cluster <-FindClusters(pbmc_Neighbor_PCA,resolution = 3 )
DotPlot(pbmc_Neighbor_PCA_cluster,features = c("Myc","Irf4","Ly75","Kdm6b","Mtor"))
