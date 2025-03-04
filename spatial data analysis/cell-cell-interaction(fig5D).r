#download known ligand-receptor pairs online
lr_network <- readRDS("lr_network_mouse_21122021.rds") #from NicheNet 
lr_network$lr_pair <- paste0(lr_network$from,'_',lr_network$to)
ligand_receptor <- readRDS('mouse_lr_pair.rds') #from CellTalk database

#calculate L-R pairs weights on spatial
library(Seurat)
seu_region <- readRDS('seu_region.rds') #interface information was preserved in 'identify_myo-macro_interface.r' step

ligand_receptor <- lr_network
lr_scale_dataset1 <- function(obj) {
  exp = GetAssayData(obj,slot='counts')
  #exp <- data.frame(matrix(exp))
  nm = rownames(exp)
  lr = ligand_receptor[apply(ligand_receptor, 1, function(x) { length(intersect(nm, x))}) == 2,]
  
  #exp = apply(exp,1,function(x) x/max(x))
  #exp = t(exp)
  exp1 = exp[lr$from,]
  exp2 = exp[lr$to,]
  exp3 = log10(exp1) + log10(exp2)
  exp3 = 10**exp3
  exp3[is.na(exp3)] = 0
  rownames(exp3) = paste(lr$from, lr$to,sep="_")
  return(exp3[rownames(exp3)[rowSums(exp3)>0],])
}
lr_scale_dataset2 <- function(obj) {
  exp = GetAssayData(obj,slot='counts')
  #exp <- data.frame(matrix(exp))
  nm = rownames(exp)
  lr = ligand_receptor[apply(ligand_receptor, 1, function(x) { length(intersect(nm, x))}) == 2,]
  
  #exp = apply(exp,1,function(x) x/max(x))
  #exp = t(exp)
  exp1 = exp[lr$ligand_gene_symbol,]
  exp2 = exp[lr$receptor_gene_symbol,]
  exp3 = log10(exp1) + log10(exp2)
  exp3 = 10**exp3
  exp3[is.na(exp3)] = 0
  rownames(exp3) = paste(lr$ligand_gene_symbol, lr$receptor_gene_symbol,sep="_")
  return(exp3[rownames(exp3)[rowSums(exp3)>0],])
}
ligand_receptor <- lr_network
seu_region_scaled1 <- lr_scale_dataset1(seu_region)
ligand_receptor <- readRDS('/data/work/spatial_riga_batch2/12_macro_sub/mouse_lr_pair.rds')
seu_region_scaled2 <- lr_scale_dataset2(seu_region)
a <- intersect(rownames(seu_region_scaled1),rownames(seu_region_scaled2))
b <- setdiff(rownames(seu_region_scaled1),rownames(seu_region_scaled2))
c <- setdiff(rownames(seu_region_scaled2),rownames(seu_region_scaled1))
seu_region_scaled1 <- seu_region_scaled1[rownames(seu_region_scaled1) %in% c(a,b),]
seu_region_scaled2 <- seu_region_scaled2[rownames(seu_region_scaled2) %in% c,]
seu_region_scaled <- rbind(seu_region_scaled1,seu_region_scaled2)

lr_seu <- CreateSeuratObject(seu_region_scaled,meta.data = seu_region@meta.data)
lr_seu <- subset(lr_seu, group == 'LN')
Idents(lr_seu) <- 'selected'
degs <- FindMarkers(lr_seu,ident.1 = 'selected',ident.2 = 'inner_medulla',only.pos=T)
degs <- degs[degs$p_val_adj < 0.05,]
degs$ligand <- sub(rownames(degs),pattern = "-.*",replacement = "")
degs$receptor <- sub(rownames(degs),pattern = ".*-",replacement = "")
degs$weights <- degs$pct.1/degs$pct.2

library(RColorBrewer)
p <- ggplot(head(degs,100),aes(x=ligand,y=receptor))+
  geom_tile(aes(fill=weights),color = "white") +
  theme_classic()+ scale_fill_gradientn(name='interaction_weight',colours=rev(brewer.pal(8, "Spectral")) ,
                           breaks=c(1,4),labels=c("Min","Max"),na.value = "transparent",
                           limits=c(1,4))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(legend.position="bottom")
p
