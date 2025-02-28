setwd('work_path')
glom_list <- list.files(pattern='cluster_label.txt') # data from 'identify_gloms.py' 
seu <- readRDS('seu_region.rds')
pods <- subset(seu,Annotation == 'Podocyte')

#integrate glom gene expression files 
for(j in c(1:6)){
a <- read.table(glom_list[j])
sample_id <- strsplit(glom_list,split='_')[[j]][1]
seu_sample <- subset(pods,sample == as.character(sample_id))
table(paste0(seu_sample$coor_x,'_',seu_sample$coor_y)==paste0(a$V1,'_',a$V2))
seu_sample$glom_label <- a$V3
    
assay <- seu_sample@assays$RNA@counts
    assay <- data.frame(assay)
    assay <- data.frame(t(assay))
    assay$glom_label <- seu_sample$glom_label
    length(colSums(assay))
    assay_all_glom <- data.frame(matrix(ncol=27901))
    colnames(assay_all_glom) <- colnames(assay)[1:27901]
    
for(i in sort(unique(seu_sample$glom_label))){
assay_glom <- assay[assay$glom_label ==i,]
assay_glom <- colSums(assay_glom)[1:27901]
assay_all_glom <- rbind(assay_all_glom,assay_glom)
}
assay_all_glom <- na.omit(assay_all_glom)
rownames(assay_all_glom) <- paste0(sample_id,'_glom_',sort(unique(seu_sample$glom_label)))
seu_glom <- CreateSeuratObject(t(assay_all_glom))
seu_glom$sample <- sample_id
if(j==1){seu_glom_allsample = seu_glom}else{seu_glom_allsample <- merge(seu_glom_allsample,seu_glom)}
}

seu_glom_allsample <- NormalizeData(seu_glom_allsample)
seu_glom_allsample <- ScaleData(seu_glom_allsample)
seu_glom_allsample <- FindVariableFeatures(seu_glom_allsample)
seu_glom_allsample <- RunPCA(seu_glom_allsample,dims=1:30)

#fig 6C
library(RColorBrewer)
cols <-  brewer.pal(11,"Spectral")
cols <- cols[c(1,2,3,9,10,11)]
sample_color <- cols 
names(sample_color) <- unique(seu_glom_allsample$sample_id)
DimPlot(seu_glom_allsample,group.by='sample_id',reduction='pca',cols=sample_color)

#fig 6D
Idents(seu_glom_allsample) <- 'group'
degs <- FindAllMarkers(seu_glom_allsample,only.pos=T)
degs <- degs[degs$p_val_adj < 0.05,]
