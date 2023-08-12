comp_normalized <- readRDS(file="../PCA/comp_normalized.rds")

Hierarchiral_clustering <- function(matrix,desired_number_of_groups) {
  
  # Dist matrix is constructed by row, so matrix is transposed here.
  hm.matrix <- as.matrix(t(matrix))
  
  # Calculate distance matrix (correlation matrix)
  y.dist <- amap::Dist(hm.matrix, method="spearman") 
  
  # Perform the clustering using the distances (build the dendrogram). 
  y.cluster <- hclust(y.dist, method="average")
  
  # Cut the dendrogram to get desired_number_of_groups
  cluster_labels <- cutree(y.cluster, k=desired_number_of_groups)
  
  return(cluster_labels)
  
}

Hierarchiral_clustering_pca_normalized_12 <- Hierarchiral_clustering(comp_normalized,12)
