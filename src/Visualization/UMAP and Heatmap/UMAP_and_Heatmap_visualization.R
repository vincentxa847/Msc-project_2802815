# umap [version 0.2.10.0]
# ggplot2 [version 3.4.1]
# ComplexHeatmap [version 2.10.0]

library(ggplot2)

comp_normalized <- readRDS(file="comp_normalized.rds")
Matrix_normalized_excluded_removeanti <- readRDS(file="Matrix_normalized_excluded_removeanti.rds")

## Produce the umap coordinates for each genes
umap_result_pca_normalized <- umap::umap(t(comp_normalized),n_neighbors=15,n_components=2)

# saveRDS(umap_result_pca_normalized,file = "umap_coordinate.rds")
umap_coordinate <- readRDS(file = "umap_coordinate.rds")

## Using UMAP to visualize clustering information 
UMAP_visualization <- function(umap_result,cluster_labels) {
  
  # Palette with 30 colours
  color_palette <- c("#8f5362","#b96570","#d37b6d","#e0a981","#ecd09c","#d4daa1","#a3c8a4","#79b4a0",
                    "#6888a5","#706d94","#FABB6E","#FC8002","#ADDB88","#369F2D","#FAC7B3","#EE4431",
                    "#1663A9","#B4B4D5","#8481BA","#B0A875",
                    "#4A3933","#92221F","#F1ab19","#055743","#d6e1d3","#E3882F","#49548A","#514549","#B28282","#B14E53")
  # Plotted data
  plot_data <- data.frame(UMAP1=umap_result$layout[, 1],UMAP2=umap_result$layout[, 2],Cluster=cluster_labels)
  
  # Create a data frame with unique cluster labels and their corresponding UMAP coordinates
  cluster_centers <- aggregate(cbind(UMAP1, UMAP2) ~ Cluster,data=plot_data,FUN=mean)
  
  
  ggp <- ggplot(plot_data) +
    geom_point(aes(x=UMAP1,y=UMAP2,color=as.factor(Cluster)),size=3) +
    scale_color_manual(values=color_palette) +
    labs(x="UMAP 1",y="UMAP 2",colour="Clusters") + 
    theme_bw()
  
  # Add cluster numbers as text labels at the cluster center
  ggp <- ggp + geom_text(data=cluster_centers,aes(x=UMAP1,y=UMAP2,label=Cluster),color="black",size=6)
  
  return(ggp)
}

# change cluster_labels for visualization of different clustering result
UMAP_visualization(umap_coordinate,kmeans_pca_normalized_8)

## Using ComplexHeatmap, this package intergrates clustering and plotting 
library("ComplexHeatmap")

# Define the color gradient
my_palette <- colorRampPalette(c("#C25160", "#5AA4AE","#806D9A","#4F794A","#06436F"))

## plot k-means clustering on normalized TPM matrix
#' ComplexHeatmap cluster row by default, but for now only need to cluster genes (column), therefore
#' using `cluster_rows = FALSE` to turn off row clustering
#' change the column_split and anno_block when visualizing different result

Heatmap(Matrix_normalized_excluded_removeanti,name="TPM",
        column_split=kmeans_pca_normalized_30,show_column_dend=FALSE,
        cluster_rows=TRUE,show_row_dend=FALSE,
        show_row_names=FALSE,
        show_column_names=FALSE,
        top_annotation=HeatmapAnnotation(foo=anno_block(gp=gpar(fill=my_palette(30)))))
