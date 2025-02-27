library(Seurat)
sn <- readRDS('sn.rds')

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

