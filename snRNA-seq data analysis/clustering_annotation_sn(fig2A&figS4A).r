library(Seurat)
sn <- readRDS('sn.rds')
sn <- subset(sn,nCount_RNA >500 & nCount_RNA < 3500 & percent.mt < 5)

#seurat pipeline to remove batch effect
sample.list <- SplitObject(control1, split.by = "sample")
sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = sample.list)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
sn.combined <- IntegrateData(anchorset = sample.anchors)
DefaultAssay(sn.combined) <- 'integrated'
sn.combined  <- ScaleData(sn.combined , verbose = FALSE)
sn.combined  <- RunPCA(sn.combined ,verbose = FALSE)
sn.combined  <- FindNeighbors(sn.combined ,dims = 1:30, verbose = FALSE)
sn.combined  <- RunUMAP(sn.combined ,dims = 1:30, verbose = FALSE)
sn.combined  <- FindClusters(sn.combined ,verbose = FALSE, resolution = 0.6)

#find top DEGs to annotate cell clusters
DefaultAssay(sn.combined) <- 'sn.combined'
degs <- FindAllMarkers(sn.combined,only.pos=T)

#annotate celltypes based on top DEGs
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in%c(0,9), 'cluster_annotation'] <- 'TAL'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(1,3,15,16), 'cluster_annotation'] <- 'PCT'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(2), 'cluster_annotation'] <- 'PCT_S3'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(4), 'cluster_annotation'] <- 'Endo'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(5), 'cluster_annotation'] <- 'CD-PC'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(6), 'cluster_annotation'] <- 'FIB'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(7), 'cluster_annotation'] <- 'DCT'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(8), 'cluster_annotation'] <- 'CNT'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(10), 'cluster_annotation'] <- 'DTL'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(11), 'cluster_annotation'] <- 'Immune'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(12), 'cluster_annotation'] <- 'ICA'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(13), 'cluster_annotation'] <- 'ICB'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(14), 'cluster_annotation'] <- 'ATL'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(17), 'cluster_annotation'] <- 'VSM'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(19), 'cluster_annotation'] <- 'Podocyte'
sn.combined@meta.data[sn.combined@meta.data$seurat_cluster %in% c(21), 'cluster_annotation'] <- 'Dividing_cell'

#fig 2A
p <- DimPlot(sn.combined,group.by='cluster_annotation',label=T)+labs(title='Annotation')+scale_color_manual(values=cols)

#fig S4A
library(ggplot2)
DefaultAssay(sn.combined) <-"RNA"
gene<-c('Slc12a1','Slc34a1','Miox','Kap','Slc13a3','Egfl7','Flt1','Fxyd4','Aqp2','Fbln5','Col1a2',
        'Slc12a3','Tmem52b','Slc8a1','Calb1','Slc4a11','Bst1','Ptprc','Slc4a9','Slc4a1','Slc26a4',
        'Cryab','S100a6','Myh11','Ren1','Nphs2','Synpo','Mki67','Top2a') #selected celltype specific genes
p<-DotPlot(sn.combined,assay="RNA", features = gene, cols = c("#E7E6E6", "red"), dot.scale = 6,col.min=0,col.max=2,dot.min=0)+coord_flip()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+coord_fixed()
p

