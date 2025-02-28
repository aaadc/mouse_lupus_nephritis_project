library(Seurat)
immune <- subset(sn.combined,seurat_clusters == 11)
immune <- CreateSeuratObject(immune@assays$RNA@counts,meta.data=immune@meta.data)

#seurat pipline to remove batch effect
sample.list <- SplitObject(immune, split.by = "sample")
sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = sample.list)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
sn.immune <- IntegrateData(anchorset = sample.anchors,k.weight=20)
DefaultAssay(sn.immune) <- 'integrated'
sn.immune  <- ScaleData(sn.immune, verbose = FALSE)
sn.immune <- RunPCA(sn.immune,verbose = FALSE)
sn.immune  <- FindNeighbors(sn.immune ,dims = 1:30, verbose = FALSE)
sn.immune  <- RunUMAP(sn.immune ,dims = 1:30, verbose = FALSE)
sn.immune  <- FindClusters(sn.immune ,verbose = FALSE)

#find top DEGs to annotate cell clusters
DefaultAssay(sn.immune) <- 'RNA'
degs <- FindAllMarkers(sn.immune,only.pos=T)

#annotate cell clusters 
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters == 0,'cluster_annotation'] = 'Tcell'
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters == 1,'cluster_annotation'] = 'Itgax_macrophage'
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters %in% c(2,3),'cluster_annotation'] = 'tubular'
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters == 4,'cluster_annotation'] = 'Bcell'
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters == 5,'cluster_annotation'] = 'Lyz2_macrophage'
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters == 6,'cluster_annotation'] = 'Endo'
sn.immune@meta.data[sn.immune@meta.data$seurat_clusters == 9,'cluster_annotation'] = 'DC'

#Fig 2B
p <- DimPlot(sn.immune,group.by='cluster_annotation',label=T)
p

#Fig 2C
library(ggplot2)
DefaultAssay(sn.immune) <-"RNA"
gene<-c('Itgax','Itgam','Lyz2','Apoe','Lck','Lat','Siglecg','Cd79b','Clec9a')
p<-DotPlot(subset(sn.immune,cluster_annotation %in% c('Tcell','Itgax_macrophage','Bcell','Lyz2_macrophage','DC')),
           assay="RNA", features = gene, cols = c("#E7E6E6", "red"), dot.scale = 6,col.min=0,col.max=2,dot.min=0)+coord_flip()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+coord_fixed()

