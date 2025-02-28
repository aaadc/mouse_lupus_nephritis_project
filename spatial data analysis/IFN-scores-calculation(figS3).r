# get IFN scores on spatial
list <- list.files(pattern='txt') #GO_term_summary* files download from database
library(Matrix)
library(data.table)
genes_all <- c()
for(i in c(1:length(list))){
gene <- data.frame(fread(list[i]))
description <- strsplit(gene$Qualifier,split=' ')
for(j in c(1:nrow(gene))){
if(c('negative') %in% unlist(description[j])){gene[j,'preserve']='true'}else{gene[j,'preserve']='false'}
}
gene_p <- gene[gene$preserve == 'false','MGI.Gene.Marker.ID']
gene_p <- unique(gene_p)
genes_all <- unique(c(genes_all,gene_p))
}

#pipeline for aucell
library(ggplot2)
library(Seurat)
library(AUCell)
seu_region <- readRDS('seu_region.rds')
exprMatrix <- seu_region@assays$RNA@counts
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
a <- list(genes_all);names(a) <- 'IFN-I'
cells_AUC <- AUCell_calcAUC(a, cells_rankings)
seu_region@meta.data <- cbind(seu_region@meta.data,t(getAUC(cells_AUC)))

p <- ggplot(seu_region@meta.data,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),fill=as.numeric(seu_region@meta.data[,'IFN-I'])))+geom_tile()+
  scale_fill_gradientn(colours =c('#0c3383','#005ea3','#0a88ba','#00c199','#f2d338','#f6b132','#f28f38','#e76124','#d91e1e'))+
  theme_bw()+coord_fixed()+facet_wrap(~sample,ncol=3)+labs(fill='IFN-I')
p
