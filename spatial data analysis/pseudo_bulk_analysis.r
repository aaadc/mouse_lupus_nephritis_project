library(Seurat)
library("DESeq2")

spatial <- readRDS('seu_region.rds')
DefaultAssay(spatial) <- 'RNA'

# get inner medulla region gene expression pattern
inner <- subset(spatial,region=='inner_medulla')
bulk_allsample <- data.frame(matrix(nrow=nrow(inner)))
for (i in unique(inner$sample)){
    sample <- subset(inner,sample == i)
    assays <- sample@assays$RNA@counts
    bulk <- rowSums(assays)
    bulk_allsample <- cbind(bulk_allsample,'counts'=bulk)
    colnames(bulk_allsample)[ncol(bulk_allsample)] <- i
    }
bulk_allsample <- bulk_allsample[,-1]
colnames(bulk_allsample) <- paste0(colnames(bulk_allsample),"_inner_medulla")
bulk_allsample_inner_medulla <- bulk_allsample
meta <- data.frame(matrix(ncol=2,nrow=6));colnames(meta) <- c('sample','group')
meta$sample = colnames(bulk_allsample)
meta$group = c(rep('LN',3),rep('CTRL',3))
rownames(meta) <- meta$sample
dds <- DESeqDataSetFromMatrix(countData=bulk_allsample, 
                              colData=meta, 
                              design=~group)
dds$group <- relevel(dds$group, ref="CTRL")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
diff_gene_deseq1 <-subset(res,padj < 0.05 & (log2FoldChange > 0.25 ))
diff_gene_deseq1 <- as.data.frame(diff_gene_deseq1)
diff_gene_deseq1 <- diff_gene_deseq1[order(diff_gene_deseq1$log2FoldChange,decreasing=T),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange < -0.25 ))
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$log2FoldChange,decreasing=F),]
diff_gene_deseq <- rbind(diff_gene_deseq1,diff_gene_deseq2)
write.csv(diff_gene_deseq1,'bulk_inner_medulla_up_imq.csv')
write.csv(diff_gene_deseq2,'bulk_inner_medulla_up_control.csv')
normalized.data <- counts(dds,normalized = TRUE)
write.csv(normalized.data,'inner_medulla_normalized_counts.csv')

# get outer medulla region gene expression pattern
outer <- subset(spatial,region=='outer_medulla')
bulk_allsample <- data.frame(matrix(nrow=nrow(outer)))
for (i in unique(outer$sample)){
    sample <- subset(outer,sample == i)
    assays <- sample@assays$RNA@counts
    bulk <- rowSums(assays)
    bulk_allsample <- cbind(bulk_allsample,'counts'=bulk)
    colnames(bulk_allsample)[ncol(bulk_allsample)] <- i
    
}
bulk_allsample <- bulk_allsample[,-1]
colnames(bulk_allsample) <- paste0(colnames(bulk_allsample),"_outer_medulla")
bulk_allsample_outer_medulla <- bulk_allsample
meta <- data.frame(matrix(ncol=2,nrow=6));colnames(meta) <- c('sample','group')
meta$sample = colnames(bulk_allsample)
meta$group = c(rep('LN',3),rep('CTRL',3))
rownames(meta) <- meta$sample
dds <- DESeqDataSetFromMatrix(countData=bulk_allsample, 
                              colData=meta, 
                              design=~group)
dds$group <- relevel(dds$group, ref="CTRL")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
diff_gene_deseq1 <-subset(res,padj < 0.05 & (log2FoldChange > 0.25 ))
diff_gene_deseq1 <- as.data.frame(diff_gene_deseq1)
diff_gene_deseq1 <- diff_gene_deseq1[order(diff_gene_deseq1$log2FoldChange,decreasing=T),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange < -0.25 ))
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$log2FoldChange,decreasing=F),]
diff_gene_deseq <- rbind(diff_gene_deseq1,diff_gene_deseq2)
write.csv(diff_gene_deseq1,'bulk_outer_medulla_up_imq.csv')
write.csv(diff_gene_deseq2,'bulk_outer_medulla_up_control.csv')
normalized.data <- counts(dds,normalized = TRUE)
write.csv(normalized.data,'outer_medulla_normalized_counts.csv')

# get cortex region gene expression pattern
cortex <- subset(spatial,region=='cortex')
bulk_allsample <- data.frame(matrix(nrow=nrow(cortex)))
for (i in unique(cortex$sample)){
    sample <- subset(cortex,sample == i)
    assays <- sample@assays$RNA@counts
    bulk <- rowSums(assays)
    bulk_allsample <- cbind(bulk_allsample,'counts'=bulk)
    colnames(bulk_allsample)[ncol(bulk_allsample)] <- i
}
bulk_allsample <- bulk_allsample[,-1]
colnames(bulk_allsample) <- paste0(colnames(bulk_allsample),"_cortex")
bulk_allsample_cortex <- bulk_allsample
meta <- data.frame(matrix(ncol=2,nrow=6));colnames(meta) <- c('sample','group')
meta$sample = colnames(bulk_allsample)
meta$group = c(rep('LN',3),rep('CTRL',3))
rownames(meta) <- meta$sample
dds <- DESeqDataSetFromMatrix(countData=bulk_allsample, 
                              colData=meta, 
                              design=~group)
dds$group <- relevel(dds$group, ref="CTRL")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
diff_gene_deseq1 <-subset(res,padj < 0.05 & (log2FoldChange > 0.25 ))
diff_gene_deseq1 <- as.data.frame(diff_gene_deseq1)
diff_gene_deseq1 <- diff_gene_deseq1[order(diff_gene_deseq1$log2FoldChange,decreasing=T),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange < -0.25 ))
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
diff_gene_deseq2 <- diff_gene_deseq2[order(diff_gene_deseq2$log2FoldChange,decreasing=F),]
diff_gene_deseq <- rbind(diff_gene_deseq1,diff_gene_deseq2)
write.csv(diff_gene_deseq1,'bulk_cortex_up_imq.csv')
write.csv(diff_gene_deseq2,'bulk_cortex_up_control.csv')
normalized.data <- counts(dds,normalized = TRUE)
write.csv(normalized.data,'cortex_normalized_counts.csv')

"fig 1E"

library(ggplot2)
library(ggrepel)
# get pca embedings 
bulk_allsample <- cbind(bulk_allsample_cortex,bulk_allsample_outer_medulla)
bulk_allsample <- cbind(bulk_allsample,bulk_allsample_inner_medulla)
meta <- data.frame(matrix(ncol=2,nrow=18));colnames(meta) <- c('sample','group')
meta$sample = substr(colnames(bulk_allsample),1,6)
meta$group = c(rep('LN',3),rep('CTRL',3))
rownames(meta) <- colnames(bulk_allsample)
meta[meta$sample == 'M09902','sample'] <- 'IMQ_1'
meta[meta$sample == 'M09903','sample'] <- 'IMQ_2'
meta[meta$sample == 'M09904','sample'] <- 'IMQ_3'
meta[meta$sample == 'M09907','sample'] <- 'CTRL_1'
meta[meta$sample == 'M09908','sample'] <- 'CTRL_2'
meta[meta$sample == 'M09909','sample'] <- 'CTRL_3'
meta[meta$sample %in% c('IMQ_1','IMQ_2','IMQ_3'),'group'] <- 'IMQ'
dds <- DESeqDataSetFromMatrix(countData=bulk_allsample, 
                              colData=meta, 
                              design=~group)
dds$group <- relevel(dds$group, ref="CTRL")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
rld <- rlog(dds, blind=TRUE)
pcadata <- plotPCA(rld, intgroup="group",returnData=TRUE,ntop=2000)
pcadata$region <- substring(pcadata$name,8)
pcadata$sample <- substring(rownames(pcadata),1,6)
pcadata$sample <- meta$sample

#calculate percentage of PC1 & PC2
rv <- rowVars(assay(rld))
ntop=2000
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
percentVa[1]*100  # Percentage of variance for PC1, 52.8 in this case
percentVar[2]*100  # Percentage of variance for PC2, 20.7 in this case

#draw the plot
p <- ggplot(pcadata,aes(x=PC1,y=PC2,color=group,shape=region))+geom_point(size=4)+theme_bw()+
       scale_color_manual(values=c("#2f5688","#CC0000"))+
       geom_text_repel(aes(label = sample), size = 4,color='black')+coord_fixed()+labs(x='PC1(52.8%)',y='PC2(20.7%)')

"fig 1F"
library(pheatmap)
library(ComplexHeatmap)

normalized <- read.csv('inner_medulla_normalized_counts.csv')
rownames(normalized) <- normalized$X
normalized <- normalized[,-1]

#get highly expressed genes in IMQ and control groups
bulk_imq_down_im <- read.csv('bulk_inner_medulla_up_control.csv') #save above
bulk_imq_down_im$region <- 'inner_medulla'
bulk_imq_down_im <- bulk_imq_down_im[order(bulk_imq_down_im$padj),]
bulk_imq_up_im <- read.csv('bulk_inner_medulla_up_imq.csv') #save above
bulk_imq_up_im$region <- 'inner_medulla'
bulk_imq_up_im <- bulk_imq_up_im[order(bulk_imq_up_im$padj),]
bulk_imq_down_im <- bulk_imq_down_im[bulk_imq_down_im$log2FoldChange < -2,]
bulk_imq_up_im <- bulk_imq_up_im[bulk_imq_up_im$log2FoldChange > 2,]
gene_inner  <- c(bulk_imq_down_im$X,bulk_imq_up_im$X)

library(stringr) # remove genes start with 'Rik' or 'Gm'
a <- gene_inner[str_ends(gene_inner,'Rik')]
b <- gene_inner[str_starts(gene_inner,'Gm')]
c <- which(gene_inner %in% a)
d <- which(gene_inner %in% b)
gene_inner <- gene_inner[-c(c,d)]
inner_p_mat_1 <- normalized[rownames(normalized) %in% bulk_imq_down_im$X,]
inner_p_mat_2 <- normalized[rownames(normalized) %in% bulk_imq_up_im$X,]
inner_p_mat <- rbind(inner_p_mat_1,inner_p_mat_2)
colnames(inner_p_mat) <- substr(colnames(inner_p_mat),1,6)
annotation_col = data.frame(Group = c(rep(c('IMQ'),3),rep(c('CTRL'),3)))
rownames(annotation_col) = colnames(inner_p_mat)
colnames(inner_p_mat) <- c('IMQ_1','IMQ_2','IMQ_3','CTRL_1','CTRL_2','CTRL_3')
rownames(annotation_col) <- colnames(inner_p_mat)

ann_colors = list(
  Group = c(CTRL="#2f5688", IMQ="#CC0000"))
genes <- c('Tgfbi','Ptprc','Itgam','Adgre1','Mpeg1','C1qa','C1qb','C1qc','Ccr5','Cx3cr1','Ccr2','Cxcl9','Cxcl13','Ikbke')

pheatmap <- pheatmap(inner_p_mat,scale='row',annotation_col=annotation_col,border=FALSE,annotation_colors = ann_colors) + 
rowAnnotation(link = anno_mark(at = which(rownames(inner_p_mat) %in% genes),
                                      labels = genes, labels_gp = gpar(fontsize = 10)))

"fig 1G"
library(clusterProfiler)
library(org.Mm.eg.db)
inner_unique_go <- enrichGO(bulk_imq_up_im$X,OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL",pvalueCutoff = 1) # you can get the input file for fig 1E step
write.csv(inner_unique_go@result,'inner_unique_go.csv')
go_pick <- head(inner_unique_go@result[inner_unique_go@result$`p.adjust` <0.05,],30)
go_pick <- go_pick[order(go_pick$p.adjust),]
go_pick$order <- -log10(go_pick$p.adjust)
go_pick <- go_pick[order(go_pick$p.adjust,decreasing = TRUE),]
go_plot <- ggplot(go_pick, aes(x=-log10(p.adjust), y=factor(Description,levels=go_pick$Description), color = Count)) +
  geom_point(aes(size=-log10(p.adjust)))+
  scale_color_gradient(low = 'blue', high = 'red', name = "gene_number")+
  labs(x= "-log(p_adj)", y="GO term")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10))+
  theme(legend.title = element_text(size = 10), 
               legend.text = element_text(size = 10))+
  theme(axis.title.y=element_blank())
