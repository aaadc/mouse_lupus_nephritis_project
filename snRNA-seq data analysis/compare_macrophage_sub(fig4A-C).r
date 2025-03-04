library(Seurat)
immune <- readRDS('../sn.immune.rds')

#find DEGs between two macrophage subclusters
DefaultAssay(immune) <- 'RNA'
degs <- FindMarkers(immune,ident.1=5,ident.2=1)

#fig 3A
library(ggplot2)
library(ggrepel)
degs$log2FC <- degs$avg_log2FC
cut_off_pvalue = 0.05
cut_off_logFC = 0.25
degs$group <- "Not Significant"
degs[degs$p_val_adj<cut_off_pvalue & degs$log2FC>cut_off_logFC,"group"] <- "Lyz2_macrophage"
degs[degs$p_val_adj<cut_off_pvalue & degs$log2FC< -cut_off_logFC,"group"] <- "Itgax_macrophage"
degs <- degs[order(degs$p_val_adj),]
degs.tail <- c('Stab1','Lyz2','Apoe','Abca1','Rasa1','Gas6','Dhx9','Igf1','Hpgds','Gas7','Cd84','Mrc1','Ms4a7','Ahnak',
            'Ctsd','Trem2','Ctsb')
degs.top  <- head(rownames(degs[degs$group == "Itgax_macrophage",]),10)
degs.label <- c(degs.top,degs.tail)
degs$label <- ""
for (j in degs.label) {degs[rownames(degs) == j,"label"] = j}
degs$logP <- -log10(degs$p_val_adj)
degs[rownames(degs) %in% c(degs.tail,degs.top),'genelabels'] <- 'TRUE'
p <- ggplot(degs, aes(x=log2FC,y=logP,color=group)) +
  geom_point() +
  scale_color_manual(values=c("#2f5688","#CC0000",'gray90'))+
  geom_text_repel(aes(x=log2FC,y=logP),label = ifelse(degs$genelabels == TRUE, as.character(degs$label),""), box.padding = unit(0.45, "lines"),hjust=1) + theme(legend.title=element_blank(),text = element_text(size=20),max.overlaps=nrow(a))+ 
  labs(x = "log2FC", y = "-log10(adj)")+ theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

fig 3C
library(clusterdffiler)
library(org.Mm.eg.db)
degs <- degs[degs$p_val_adj < 0.05 & degs$avg_log2FC >0,]
lyz2_go <- enrichGO(rownames(degs),OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL",pvalueCutoff = 1)
lyz2_go <- lyz2_go@result[lyz2_go@result$`p.adjust` <0.05,]
go_pick <- lyz2_go[lyz2_go$Description %in% c('positive regulation of cholesterol efflux',
    'positive regulation of fibroblast dfliferation',
    'locomotory behavior',
    'negative regulation of interleukin-1 dfduction',
    'regulation of metal ion transport'),]
p <- ggplot(go_pick, aes(x=-log10(p.adjust), y=factor(Description,levels=go_pick$Description), color = Count)) +
  geom_point(aes(size=-log10(p.adjust)))+
  scale_color_gradient(low = 'blue', high = 'red', name = "gene_number")+
  labs(x= "-log(p_adj)", y="GO term")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 9))+
  theme(legend.key.width = unit(0.3, "cm"),legend.key.height = unit(0.3, "cm"))+
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))+
  theme(axis.title.y=element_blank())

fig 3B
# preserve homologous genes
human <- readRDS('human_ln_immune.rds') #download online dataset
human_macro <- subset(human,cluster %in% c('CM0','CM1','CM2','CM3','CM4'))
mouse <- readRDS('sn.immune.rds')
mouse <- subset(mouse,cluster_annotation %in% c('Lyz2_macrophage','Itgax_macrophage'))
assay <- mouse@assays$RNA@counts
rownames(assay) <- toupper(rownames(assay))
assay <- assay[rownames(assay) %in% rownames(human),]
assay <- data.frame(assay)
assay_human <- human_macro@assays$RNA@counts
assay_human <- assay_human[rownames(assay_human) %in% rownames(assay),]
mouse_in <- CreateSeuratObject(assay,meta.data=mouse@meta.data)
human_macro <- CreateSeuratObject(assay_human,meta.data=human_macro@meta.data)
mouse_in$source <- 'mouse'
mouse_in$cluster <- mouse_in$cluster_annotation
human_macro$source <- 'human'

# seurat label transfer pipeline
human_macro <- NormalizeData(human_macro)
human_macro <- FindVariableFeatures(human_macro)
human_macro <- ScaleData(human_macro)
human_macro <- RunPCA(human_macro)
human_macro <- FindNeighbors(human_macro, dims = 1:30)
human_macro <- FindClusters(human_macro)
mouse_in <- NormalizeData(mouse_in)
mouse_in <- FindVariableFeatures(mouse_in)
mouse_in <- ScaleData(mouse_in)
mouse_in <- RunPCA(mouse_in)
mouse_in <- FindNeighbors(mouse_in, dims = 1:30)
mouse_in <- FindClusters(mouse_in)
anchors <- FindTransferAnchors(reference = human_macro, query = mouse_in, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = human_macro$cluster, dims = 1:30)

mouse_in <- AddMetaData(mouse_in, metadata = predictions)
library(dplyr)
df <- mouse_in@meta.data %>% group_by(cluster, predicted.id) %>% summarize('sum'= sum(prediction.score.max))

df[df$predicted.id == 'CM0','predicted'] <- 'CM0 (Inflammatory CD16+ macrophages)'
df[df$predicted.id == 'CM1','predicted'] <- 'CM1 (Phagocytic CD16+ macrophages)'
df[df$predicted.id == 'CM2','predicted'] <- 'CM2 (Tissue resident macrophages)'
df[df$predicted.id == 'CM3','predicted'] <- 'CM3 (cDCs)'
df[df$predicted.id == 'CM4','predicted'] <- 'CM4 (M2-like CD16+ macrophages)'

p <- ggplot(df,aes(x=cluster,y=sum,fill=predicted))+geom_bar(position="fill", stat="identity")+
theme_classic()+labs(x=NULL,y='proportion',fill=NULL)+
theme(axis.text.x = element_text(size=16, vjust = 0.5))+
theme(axis.text.y = element_text(size=16, vjust = 0.5))+
theme(legend.position = "bottom")+
  guides(fill = guide_legend(ncol = 4))+ 
scale_fill_manual(values=c('#6DBB86','#B8DFB9','#76BDE5','#E7E6D6'))
p
