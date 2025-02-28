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
seu_glom_allsample@meta.data[seu_glom_allsample@meta.data$sample %in% c('M09902','M09903','M09904'),'group'] <- 'IMQ'
seu_glom_allsample@meta.data[seu_glom_allsample@meta.data$sample %in% c('M09907','M09908','M09909'),'group'] <- 'CTRL'
library(RColorBrewer)
cols <-  brewer.pal(11,"Spectral")
cols <- cols[c(1,2,3,9,10,11)]
sample_color <- cols 
names(sample_color) <- unique(seu_glom_allsample$sample_id)
DimPlot(seu_glom_allsample,group.by='sample_id',reduction='pca',cols=sample_color)

#fig 6D
library(ggrepel)
library(ggplot2)
Idents(seu_glom_allsample) <- 'group'
degs <- FindAllMarkers(seu_glom_allsample,only.pos=T)
library(stringr)
a <- degs
b <- str_ends(a$gene,'Rik')
a <- a[!a$gene %in% a[b,'gene'],]
c <- str_starts(a$gene,'Gm')
a <- a[!a$gene %in% a[c,'gene'],]
a$gene <- sub(a$gene,pattern='\\.',replacement = '-')
a$log2FC <- a$avg_log2FC
a[a$cluster == 'CTRL','log2FC'] <-  -a[a$cluster == 'CTRL','log2FC']
cut_off_pvalue = 0.05
cut_off_logFC = 0.25
a$group <- "Not Significant"
a[a$p_val_adj<cut_off_pvalue & a$log2FC>cut_off_logFC,"group"] <- "IMQ"
a[a$p_val_adj<cut_off_pvalue & a$log2FC< -cut_off_logFC,"group"] <- "CTRL"
a$logP <- -log10(a$p_val_adj)
a <- a[order(a$logP),]
#a.tail <- tail(a[a$group == "up","X"],30)
#a.top  <- head(a[a$group == "down","X"],30)
#a.label <- c(a.top,a.tail)
a.tail <- tail(a[a$group == "IMQ","gene"],20)
a.top  <- tail(a[a$group == "CTRL","gene"],20)
a.label <- c(a.top,a.tail)
a$label <- ""
for (j in a.label) {a[a$gene == j,"label"] = j}
p1 <- ggplot(data=a, aes(x=log2FC, y=logP, col=group, label=label)) + 
        geom_point() + 
        theme_bw() +
        geom_text_repel() +
        scale_color_manual(values=c("#2f5688","#CC0000",'gray')) +
        #geom_vline(xintercept=c(-0.25, 0.25), col="gray30",linetype='dashed')+
    theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())+labs(fill='Group')+theme(legend.position="none")
p1
