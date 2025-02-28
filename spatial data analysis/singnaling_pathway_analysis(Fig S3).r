#progeny pileline

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

seu_list <- list.files(pattern = '*.rds') #load spatial data

for (i in c(1:6)){
sample <- strsplit(seu_list,"_")[[i]][1]
spatial <- readRDS(seu_list[i])
spatial$sample <- sample
DefaultAssay(spatial) <- 'RNA'

Idents(spatial) <- "sample"
CellsClusters <- data.frame(Cell = names(Idents(spatial)), 
                            CellType = as.character(Idents(spatial)),
                            stringsAsFactors = FALSE)
spatial <- progeny(spatial, scale=FALSE, organism="Mouse", top=500, perm=1, return_assay = TRUE)
spatial <- Seurat::ScaleData(spatial, assay = "progeny") 
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(spatial, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

progeny_scores_df$coor_x <- sub(progeny_scores_df$Cell,pattern = "_.*",replacement = "")
progeny_scores_df$coor_x <- as.numeric(progeny_scores_df$coor_x)
progeny_scores_df$coor_x <- progeny_scores_df$coor_x - min(progeny_scores_df$coor_x)
progeny_scores_df$coor_y <- sub(progeny_scores_df$Cell,pattern = ".*_",replacement = "")
progeny_scores_df$coor_y <- as.numeric(progeny_scores_df$coor_y)
progeny_scores_df$coor_y <- progeny_scores_df$coor_y - min(progeny_scores_df$coor_y)
progeny_scores_df$sample <- sample
if(i==1){progeny_scores_df_all_sample <- progeny_scores_df}else{progeny_scores_df_all_sample <- rbind(progeny_scores_df_all_sample,progeny_scores_df)}  
pathways <- unique(progeny_scores_df$Pathway)
}

meta <- progeny_scores_df_all_sample[progeny_scores_df_all_sample$Pathway == 'NFkB',] #selected signaling pathways
p <- ggplot(meta,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),fill=Activity))+geom_tile()+
  scale_fill_gradientn(labels=c("Min","Max"),breaks=c(min(meta$Activity),max(meta$Activity)),
                       colours =c('#0c3383','#005ea3','#0a88ba','#00c199','#f2d338','#f6b132','#f28f38','#e76124','#d91e1e'))+
  theme_bw()+coord_fixed()+facet_wrap(~sample_id,ncol=6)+labs(x='x',y='y',fill='NFkB')+
    theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())+theme(
  strip.text = element_text(size = 15))
p
