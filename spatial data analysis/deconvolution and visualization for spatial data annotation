# pileline for spacexr

library("spacexr")
library("Seurat")
snRNA <- readRDS('snRNA_all_library.rds')
DefaultAssay(snRNA) <- "RNA"
ref_counts <- GetAssayData(snRNA, slot = "counts")
ref_celltypes <- snRNA@meta.data$cluster_annotation
ref_celltypes <- as.factor(ref_celltypes)
names(ref_celltypes) <- rownames(snRNA@meta.data)
nUMI <- snRNA@meta.data$nCount_RNA; names(nUMI) <- rownames(snRNA@meta.data)
reference <- Reference(ref_counts, ref_celltypes, nUMI)

spatial <- readRDS(sample_id_spatial.rds) 
spatial_counts <- GetAssayData(spatial, slot = "counts")
spatial_coords <- spatial@meta.data[,c("coor_x","coor_y")]
nUMI <- colSums(spatial_counts)
puck <- SpatialRNA(spatial_coords, spatial_counts, nUMI)
barcodes <- colnames(puck@counts)
myRCTD <- create.RCTD(puck, reference, max_cores = 1,CELL_MIN_INSTANCE = 5,UMI_min = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]]
spatialRNA <- myRCTD@spatialRNA
setwd('/your_work_path')
saveRDS(myRCTD, paste0("./", "/",'sample_id',"my_RCTD", ".rds"))
#got all RCTD annotation files for all sample this step

#fig 2D
setwd('/your_work_path')
rctd_list <- list.files(pattern = '*my_RCTD.rds')
seu_list <- list.files(pattern = '*sample_id_spatial.rds')
library(ggplot2)
library(data.table)

scores_allsample <- data.frame(matrix(ncol=1))
for(i in c(1:6)){
sample <- strsplit(seu_list,split = '_')[[i]][1] #get sampel name
myRCTD <- readRDS(rctd_list[i])
result <- myRCTD@results$results_df
spatial <- readRDS(paste0('/data/work/spatial_riga_batch2/01_02_bin50_seuobj/',seu_list[i]))
spatial$cell_id <- colnames(spatial)
spatial <- subset(spatial,cell_id %in% rownames(result))
spatial$coor_x <- spatial$coor_x - min(spatial$coor_x)
spatial$coor_y <- spatial$coor_y - min(spatial$coor_y)
scores <- as.data.frame(as.matrix(myRCTD@results$weights))
dataa <- apply(scores,1,function(x){names(scores)[which.max(x)]})
scores <- cbind(scores,spatial@meta.data[,c('coor_x','coor_y')])
scores$sample <- sample
scores$celltype <- dataa
scores_allsample <- bind_rows(scores_allsample,scores)
}
scores_allsample <- scores_allsample[-1,]

library(ggplot2)
library(ggsci)
cols <- unique(c(pal_locuszoom("default")(7),pal_igv("default")(51)))
my_cols <- cols[1:20]
names(my_cols) <- sort(unique(scores_allsample$celltype));my_cols[which(names(my_cols)=='TAL')] <- 'gray94'
p <- ggplot(scores_allsample,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),fill=celltype))+
    facet_wrap(~factor(sample),ncol=3)+geom_tile()+coord_fixed()+
    scale_fill_manual(values=my_cols)+
    theme_bw()+labs(x='x',y='y',fill='Annotation')+
    theme(legend.position="bottom")+
    theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())


#fig S5 A-E
library(ggplot2)
for(i in celltype){
p1 <- ggplot(scores,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),fill=as.numeric(scores[,i])))+geom_tile()+coord_fixed()+
    facet_wrap(~sample_id,ncol=6)+
    theme_bw()+labs(x='x',y='y')+
    scale_fill_gradient(low = "gray94", high = "red",limits = c(0, 0.3))+labs(title=i,fill='RCTD_score')+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())+theme(
    strip.text = element_text(size = 15))+
    theme(plot.title = element_text(size = 15))
print(p1)
}
  
  
