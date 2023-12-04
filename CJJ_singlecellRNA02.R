library(patchwork)
library(Seurat)
library(dplyr)
library(future)

pbmc<-readRDS("pbmc.rds")
pbmc_Totalgene<-as.data.frame(t(as.data.frame(pbmc@assays$RNA@data)[c("Myc","Irf4","Ly75","Kdm6b","Mtor"),]))
pbmc_Myc<-pbmc_Totalgene[which(pbmc_Totalgene$Myc > 0 ),]
pbmc_Irf4<-pbmc_Totalgene[which(pbmc_Totalgene$Irf4 > 0),]
pbmc_Ly75<-pbmc_Totalgene[which(pbmc_Totalgene$Ly75 > 0),]
pbmc_Kdm6b<-pbmc_Totalgene[which(pbmc_Totalgene$Kdm6b > 0),]
pbmc_Mtor<-pbmc_Totalgene[which(pbmc_Totalgene$Mtor > 0),]

write.table(pbmc_Irf4,"pbmc_Irf4.txt",sep = "\t")
write.table(pbmc_Kdm6b,"pbmc_Kdm6b.txt",sep = "\t")
write.table(pbmc_Ly75,"pbmc_Ly75.txt",sep = "\t")
write.table(pbmc_Mtor,"pbmc_Mtor.txt",sep = "\t")
write.table(pbmc_Myc,"pbmc_Myc.txt",sep = "\t")

remove(pbmc_Totalgene)
pbmc_Irf4<-rownames(pbmc_Irf4)
pbmc_Kdm6b<-rownames(pbmc_Kdm6b)
pbmc_Ly75<-rownames(pbmc_Ly75)
pbmc_Mtor<-rownames(pbmc_Mtor)
pbmc_Myc<-rownames(pbmc_Myc)

pbmc@meta.data$Myc <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Irf4 <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Ly75 <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Kdm6b <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Mtor <- rep("neg",length(rownames(pbmc@meta.data)))

for (i in pbmc_Irf4) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Irf4[j] <- "Irf4"
    }
  }
}

for (i in pbmc_Kdm6b) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Kdm6b[j] <- "Kdm6b"
    }
  }
}

for (i in pbmc_Ly75) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Ly75[j] <- "Ly75"
    }
  }
}

for (i in pbmc_Mtor) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Mtor[j] <- "Mtor"
    }
  }
}

for (i in pbmc_Myc) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Myc[j] <- "Myc"
    }
  }
}

pbmc@meta.data$Mix<- paste(pbmc@meta.data$Myc,pbmc@meta.data$Irf4,pbmc@meta.data$Kdm6b,pbmc@meta.data$Ly75,pbmc@meta.data$Mtor,sep = "_")
pbmc_lda<-RunLDA(pbmc,labels = pbmc$Mix)
pbmc_lda<-RunUMAP(pbmc_lda,reduction = "lda",reduction.name = "lda_umap",dims = 1:31)
pbmc_lda<-RunTSNE(pbmc_lda,reduction = "lda",reduction.name = "lda_tsne",dims = 1:31)
Idents(pbmc_lda)<-pbmc_lda$Mix
DotPlot(pbmc_lda,features = c("rna_Myc","rna_Irf4","rna_Kdm6b","rna_Ly75","rna_Mtor"))
library(ggplot2)
VlnPlot(pbmc_lda,features = c("rna_Myc"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda,features = c("rna_Irf4"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda,features = c("rna_Kdm6b"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda,features = c("rna_Ly75"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda,features = c("rna_Mtor"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
pbmc_lda.markers<-FindAllMarkers(pbmc_lda,only.pos = TRUE,min.pct = 0)
write.table(pbmc_lda.markers,"pbmc_lda.markers.txt",sep = "\t")                                                                                                
cellratio<-as.data.frame(proportions(table(Idents(pbmc_lda))))
cellnumber<-as.data.frame(table(Idents(pbmc_lda)))
write.table(cbind(cellratio,cellnumber),"cellratio&number.txt",sep = "\t")
saveRDS(pbmc_lda,"pbmc_lda.rds")
remove(i,j,pbmc_Irf4,pbmc_Kdm6b,pbmc_Ly75,pbmc_Mtor,pbmc_Myc,pbmc_lda.markers,cellnumber,cellratio)

pbmc@meta.data$Mix_sub<-paste(pbmc@meta.data$Myc,pbmc@meta.data$Irf4,pbmc@meta.data$Ly75,sep = "_")
pbmc_lda_sub<-RunLDA(pbmc,labels = pbmc$Mix_sub)
pbmc_lda_sub<-RunUMAP(pbmc_lda_sub,reduction = "lda",reduction.name = "lda_umap",dims = 1:7)
pbmc_lda_sub<-RunTSNE(pbmc_lda_sub,reduction = "lda",reduction.name = "lda_tsne",dims = 1:7)
Idents(pbmc_lda_sub)<-pbmc_lda_sub$Mix_sub
DotPlot(pbmc_lda_sub,features = c("rna_Myc","rna_Irf4","rna_Kdm6b","rna_Ly75","rna_Mtor"))
library(ggplot2)
VlnPlot(pbmc_lda_sub,features = c("rna_Myc"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda_sub,features = c("rna_Irf4"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda_sub,features = c("rna_Kdm6b"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda_sub,features = c("rna_Ly75"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
VlnPlot(pbmc_lda_sub,features = c("rna_Mtor"),pt.size = 0,sort = TRUE)  + theme(axis.text.x = element_text(angle = 90))
pbmc_lda_sub.markers<-FindAllMarkers(pbmc_lda_sub,only.pos = TRUE,min.pct = 0)
write.table(pbmc_lda_sub.markers,"pbmc_lda_sub.markers.txt",sep = "\t")                                                                                                
cellratio<-as.data.frame(proportions(table(Idents(pbmc_lda_sub))))
cellnumber<-as.data.frame(table(Idents(pbmc_lda_sub)))
write.table(cbind(cellratio,cellnumber),"cellratio&number.txt",sep = "\t")
saveRDS(pbmc_lda_sub,"pbmc_lda_sub.rds")

Myc<-FetchData(pbmc_lda,vars = "rna_Myc")
idents<-Idents(pbmc_lda)[colnames(pbmc_lda)]
Myc$ident<-idents
noise<-rnorm(n=length(x=Myc[,"rna_Myc"]))/1e+05
Myc[,"rna_Myc"]<-Myc[,"rna_Myc"]+noise
Myc$ident<-factor(Myc$ident,levels = names(x=rev(x=sort(x=tapply(X=Myc[,"rna_Myc"],INDEX = Myc$ident,FUN = mean),decreasing = FALSE))))

Irf4<-FetchData(pbmc_lda,vars = "rna_Irf4")
idents<-Idents(pbmc_lda)[colnames(pbmc_lda)]
Irf4$ident<-idents
noise<-rnorm(n=length(x=Irf4[,"rna_Irf4"]))/1e+05
Irf4[,"rna_Irf4"]<-Irf4[,"rna_Irf4"]+noise
Irf4$ident<-factor(Irf4$ident,levels = names(x=rev(x=sort(x=tapply(X=Irf4[,"rna_Irf4"],INDEX = Irf4$ident,FUN = mean),decreasing = FALSE))))

Ly75<-FetchData(pbmc_lda,vars = "rna_Ly75")
idents<-Idents(pbmc_lda)[colnames(pbmc_lda)]
Ly75$ident<-idents
noise<-rnorm(n=length(x=Ly75[,"rna_Ly75"]))/1e+05
Ly75[,"rna_Ly75"]<-Ly75[,"rna_Ly75"]+noise
Ly75$ident<-factor(Ly75$ident,levels = names(x=rev(x=sort(x=tapply(X=Ly75[,"rna_Ly75"],INDEX = Ly75$ident,FUN = mean),decreasing = FALSE))))

Kdm6b<-FetchData(pbmc_lda,vars = "rna_Kdm6b")
idents<-Idents(pbmc_lda)[colnames(pbmc_lda)]
Kdm6b$ident<-idents
noise<-rnorm(n=length(x=Kdm6b[,"rna_Kdm6b"]))/1e+05
Kdm6b[,"rna_Kdm6b"]<-Kdm6b[,"rna_Kdm6b"]+noise
Kdm6b$ident<-factor(Kdm6b$ident,levels = names(x=rev(x=sort(x=tapply(X=Kdm6b[,"rna_Kdm6b"],INDEX = Kdm6b$ident,FUN = mean),decreasing = FALSE))))

data<-data.frame(levels(Myc$ident),levels(Irf4$ident),levels(Ly75$ident),levels(Kdm6b$ident))
write.table(data,"data.txt",sep = "\t")
