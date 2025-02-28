library(Seurat)
scores_allsample <- read.csv('rctd_scores_allsample.csv')
rownames(scores_allsample) <- scores_allsample[,1]
scores_allsample <- scores_allsample[,-1]
seu_region <- readRDS('seu_region.rds')
seu_region$Lyz2_macrophage <- scores$Lyz2_macrophage
seu_region$FIB <- scores$FIB

#identify myofibroblasts
seu_region$Acta2 <- seu_region@assays$RNA@counts['Acta2',]
seu_region$cor1 <-  seu_region@assays$RNA@counts['Acta2',] >0 & seu_region@meta.data$FIB > 0.1
seu_region$cor2 <-  seu_region$Lyz2_macrophage > 0.1
#seu_region@meta.data[seu_region@meta.data$sample == 'M09904','cor2'] <-  seu_region@meta.data[seu_region@meta.data$sample == 'M09904','Lyz2_macrophage'] > 0.15
seu_region$cor3 <- seu_region$cor1 == 'TRUE' & seu_region$cor2 == 'TRUE'
seu_region$cor <- 'other'
seu_region@meta.data[seu_region@meta.data$cor1 == 'TRUE','cor'] <- 'myopfibrblast'
seu_region@meta.data[seu_region@meta.data$cor2 == 'TRUE','cor'] <- 'Lyz2_macrophage'
seu_region@meta.data[seu_region@meta.data$cor2 == 'TRUE'& seu_region@meta.data$cor1 == 'TRUE','cor'] <- 'co-expression'
seu_region$cor4 <-  seu_region$TAL >0.1
meta <- seu_region@meta.data

#fig 3D
p1 <- ggplot(meta, aes(as.numeric(coor_x),as.numeric(coor_y),fill = cor4))+                                        
  geom_tile()+
  theme_bw()+
  coord_fixed()+
  labs(x='x',y='y')+scale_fill_manual(values=c('gray96','gray'))+
  theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  facet_wrap(~sample,ncol=6)+scale_x_continuous(limits = c(0, 205))+
  scale_y_continuous(limits = c(0, 202))

p2 <- p1 + geom_tile(data =meta[meta$cor2 == 'TRUE',] ,aes(x = as.numeric(coor_x), y = as.numeric(coor_y)), fill = "blue")+
      geom_rect(xmin = 55, xmax = 85,   ymin = 75, ymax = 100,fill = NA, color = "black",size=0.5)+
      theme(legend.position="bottom")

#fig 5A
#right
p3 <- p1 + geom_tile(data =meta[meta$cor1 == 'TRUE',] ,aes(x = as.numeric(coor_x), y = as.numeric(coor_y)), fill = "red")+
theme(legend.position="bottom")
#left
p4 <- p2 + 
      geom_tile(data =meta[meta$cor1 == 'TRUE',] ,aes(x = as.numeric(coor_x), y = as.numeric(coor_y)), fill = "red")+
      geom_tile(data =meta[meta$cor2 == 'TRUE',] ,aes(x = as.numeric(coor_x), y = as.numeric(coor_y)), fill = "blue")+
      geom_tile(data =meta[meta$cor3 == 'TRUE',] ,aes(x = as.numeric(coor_x), y = as.numeric(coor_y)), fill = "green")+
      geom_rect(xmin = 55, xmax = 85,   ymin = 75, ymax = 100,fill = NA, color = "black",size=0.5)+
      theme(legend.position="bottom")+
      theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#fig 5B
library(RANN)
neighbor <- function(df,sample_id){
meta <- df[df$sample_id == sample_id,]
b <- meta[meta$cor1 == 'TRUE'&meta$region == 'inner_medulla',c('coor_x','coor_y')]
c <- meta[meta$cor2 == 'TRUE'&meta$region == 'inner_medulla',c('coor_x','coor_y')]
nearest <- nn2(b, c,k=1)  
nearest.1 <- data.frame(matrix(unlist(nearest), ncol=length(nearest))) 
nearest.2 <- data.frame(table(nearest.1$X2))
d <- b[nearest.1[nearest.1$X2 <= sqrt(2),'X1'],]
nearest <- nn2(c, b,k=1)  
nearest.3 <- data.frame(matrix(unlist(nearest), ncol=length(nearest))) 
nearest.4 <- data.frame(table(nearest.3$X2))
e <- c[nearest.3[nearest.3$X2 <= sqrt(2),'X1'],]
neighbor <- nearest.4
return(neighbor)
}
imq1_neighbor <- neighbor(seu_region@meta.data,'IMQ_1')
imq2_neighbor <- neighbor(seu_region@meta.data,'IMQ_2')
imq3_neighbor <- neighbor(seu_region@meta.data,'IMQ_3')

pie <- merge(imq1_neighbor,imq2_neighbor,by='Var1',all = TRUE)
pie <- merge(pie,imq3_neighbor,by='Var1',all = TRUE)
pie[is.na(pie)] <- 0
pie_plot <- data.frame(matrix(nrow=3,ncol=1))
colnames(pie_plot) <- 'number'
rownames(pie_plot) <- c('co-exist','neighbor','other')
sum <- rowSums(pie[,2:4])
pie_plot[1,1] <- sum[1]
pie_plot[2,1] <- sum(sum[2:3])
pie_plot[3,1] <- sum(sum[4:19])
pie_plot$state <- rownames(pie_plot)
pie_plot$percentage <- pie_plot$number/sum(pie_plot$number)
p1 <- ggplot(pie_plot, aes(x="", y=number, fill=state)) +
  geom_bar(stat="identity", width=1, color="white",alpha=0.8) +
  coord_polar("y", start=0) + scale_fill_brewer(palette="Set3")+theme_void()
p <- p1 + geom_text(
    aes(label = percentage),
    position = position_stack(vjust = 0.5),color='black',size=6
  )

fig 5C
library(clusterProfiler)
library(org.Mm.eg.db)

neighbor <- function(df,sample_id){
meta <- df[df$sample == sample_id,]
b <- meta[meta$cor1 == 'TRUE'&meta$region == 'inner_medulla',c('coor_x','coor_y')]
c <- meta[meta$cor2 == 'TRUE'&meta$region == 'inner_medulla',c('coor_x','coor_y')]
nearest <- nn2(b, c,k=1)  
nearest.1 <- data.frame(matrix(unlist(nearest), ncol=length(nearest))) 
nearest.2 <- data.frame(table(nearest.1$X2))
d <- b[nearest.1[nearest.1$X2 <= sqrt(2),'X1'],]
nearest <- nn2(c, b,k=1)  
nearest.3 <- data.frame(matrix(unlist(nearest), ncol=length(nearest))) 
nearest.4 <- data.frame(table(nearest.3$X2))
e <- c[nearest.3[nearest.3$X2 <= sqrt(2),'X1'],]
neighbor <- nearest.4
return(c(rownames(d),rownames(e)))
}
library(RANN)
imq1_neighbor <- neighbor(seu_region@meta.data,'M09902')
imq2_neighbor <- neighbor(seu_region@meta.data,'M09903')
imq3_neighbor <- neighbor(seu_region@meta.data,'M09904')
seu_region@meta.data[rownames(seu_region@meta.data) %in% c(imq1_neighbor,imq2_neighbor,imq3_neighbor),'preserve'] <- 'TRUE'
seu_region@meta.data$selected <- seu_region@meta.data$region
seu_region@meta.data[seu_region@meta.data$preserve %in% 'TRUE','selected'] <- 'selected'

Idents(seu_region) <- 'selected'
degs <- FindMarkers(seu_region,ident.1='selected',ident.2='inner_medulla',only.pos=T)
degs <- degs[degs$p_val_adj < 0.05,]
go <- enrichGO(rownames(degs),OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL",pvalueCutoff = 1)
write.csv(go@result[go@result$p.adjust < 0.05,],'lyz2-tal-interface-go-term.csv')
go <- go@result[go@result$p.adjust < 0.05,]

go_pick <- c('actin filament organization','leukocyte migration','regulation of apoptotic signaling pathway',
'epithelial cell migration','leukocyte chemotaxis','response to transforming growth factor beta',
'myeloid cell homeostasis','smooth muscle cell proliferation','extracellular matrix organization',
'mesenchyme development','mesenchymal cell differentiation','epithelial to mesenchymal transition',
'fibroblast migration','fibroblast proliferation')

go_pick <- go[go$Description %in% go_pick,]
p <- ggplot(go_pick, aes(x=-log10(p.adjust), y=factor(Description,levels=rev(go_pick$Description)), fill = Count)) +geom_bar(stat = "identity")+
  scale_fill_gradient(low = 'blue', high = 'red', name = "gene_number")+
  labs(x= "-log(p_adj)", y="GO term")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))+
  theme(axis.title.x = element_text(size=8),
  axis.title.y = element_text(size=8))
p

fig S7A
library(pheatmap)
gene_pick <- c('Flna','Col1a1','Tgfbr2','Rtn4','Tgfbi','Pdgfra','Anxa2','Lgals3','Aebp1','Ier3ip1','Fn1') #genes got from GO terms in fig5C
seu_region@meta.data[seu_region@meta.data$preserve %in% 'TRUE','selected'] <- 'Interface'
seu_region@meta.data[seu_region@meta.data$selected %in% c('inner_medulla'),'selected'] <- 'Inner_medulla'
seu_region@meta.data[seu_region@meta.data$selected %in% c('cortex','outer_medulla'),'selected'] <- 'Cortex & Outer_medulla'
mean <- AverageExpression(seu_region,assay='RNA',gene_pick)
mean <- mean$RNA
p <- pheatmap(mean,scale='row',cluster_cols=FALSE)


