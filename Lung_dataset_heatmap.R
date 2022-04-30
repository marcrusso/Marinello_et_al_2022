library(ComplexHeatmap)
rownames(Metadati_lung)<-Metadati_lung$samples
Metadati_lung<-as.data.frame(Metadati_lung[2])
LUNG_matrix_scale <- t(apply(tot_matrix_FPKMUQ_tumor,1,function(x){scale(x, center = T, scale = T)}))
colnames(LUNG_matrix_scale)<- colnames(tot_matrix_FPKMUQ_tumor)
LUNG_matrix_scale<-na.omit(LUNG_matrix_scale)
ha<-HeatmapAnnotation(df = Metadati_lung_ok)
heatmap<-Heatmap(LUNG_matrix_scale,heatmap_height = unit(20,"cm"), heatmap_width = unit(20,"cm"),,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,
        top_annotation =ha,column_split = Metadati_lung_ok)
draw(heatmap)
