library(patchwork)
library(Seurat)
library(dplyr)
library(future)

pbmc_lda<-readRDS("d:/20231201/pbmc_lda.rds")

pbmc_lda<-RunPCA(object = pbmc_lda)
pbmc_lda<-FindNeighbors(object = pbmc_lda,dims = 1:50)
pbmc_lda<-FindClusters(object = pbmc_lda,fea, resolution = 3)
DotPlot(pbmc_lda,features = c("rna_Myc","rna_Irf4","rna_Ly75","rna_Kdm6b"))
