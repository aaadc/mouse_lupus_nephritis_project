library(Seurat)


# find highly expresssed genes in fibroblasts in IMQ group for spatial data
spatial <- readRDS('seu_region.rds')
fib <- subset(spatial, Annotation == 'FIB' & region == 'inner_medulla')
fib <- NormalizeData(fib)
fib <- ScaleData(fib)
fib <- FindVariableFeatures(fib)
Idents(fib) <- 'group'
degs <- FindMarkers(fib,ident.1='LN',logfc.threshold =0)
spatial_fib_ln_degs <- rownames(degs[degs$avg_log2FC >0 & degs$p_val_adj <0.05,])

# find fibroblasts specific markers in snRNA-seq data
snrna <- readRDS('sn.combined.rds')
snrna_ln <- subset(snrna, group == 'imq')
DefaultAssay(snrna_ln) <- 'RNA'
snrna_ln <- NormalizeData(snrna_ln)
snrna_ln <- ScaleData(snrna_ln)
snrna_ln <- FindVariableFeatures(snrna_ln)
Idents(snrna_ln) <- 'seurat_clusters'
degs_snrna_fib <- FindMarkers(snrna_ln,ident.1=6,only.pos=T)
degs_genes <- head(rownames(degs_snrna_fib[degs_snrna_fib$p_val_adj ==0,]),200)

overlap <- intersect(spatial_fib_ln_degs,degs_genes)

#fig 3A
library(ggplot2)
library(ggrepel)
average_expression <- AverageExpression(fib,assays='RNA',rownames(degs))
average_expression <- data.frame(average_expression$RNA)
average_expression <- log10(average_expression)
average_expression$avg_log2FC <- degs$avg_log2FC
average_expression$p_val_adj <- degs$p_val_adj
average_expression$fill <- 'black'
average_expression[average_expression$p_val_adj < 0.05,'fill'] <- 'red'
for(i in overlap){average_expression[rownames(average_expression) == i,'label'] <- i}
p <- ggplot(average_expression,aes(x=LN,y=avg_log2FC,color=fill))+
geom_point(size=0.5)+scale_color_manual(values=c('gray','red'))+
geom_text_repel(aes(label=label), size=4,color='black',check_overlap = T,max.overlaps = 50)+theme_bw()+
labs(x='log2FC',y='log10(mean_expression')+
    theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())+theme(legend.position = "none")

#fig 3B
Idents(fib) <- 'group'
fib@active.ident <- factor(fib@active.ident,levels=c('CTRL','IMQ'))
p <- VlnPlot(fib,c('Dcn','Serping1','Axl','Zeb2'),cols=c("#2f5688","#CC0000"),ncol=2)
