#### k means clustering ####
# kmeans from base package stats
# cluster [version 2.1.2]
# factoextra [version 1.0.7]
# amap [version 0.8-19]

comp_normalized <- readRDS(file="../PCA/comp_normalized.rds")

## Determine the number of clusters (Elbow method) 
kmeans_Elbow <- function(data_used,number_of_clusters,nstart) {
  
  # set a empty vector to store wcss
  wcss <- vector()
  
  # set the maximun number of clusters (k) to test 
  max_k <- number_of_clusters
  
  # Run K mean clustering for different values of k and calculate WCSS
  for (k in 1:max_k) {
    
    # Set the number of runs for each k
    nstart <- nstart
    
    # Run k-means clustering
    kmeans.results <- kmeans(t(data_used), centers=k, nstart=nstart)
    
    # Extract total WCSS 
    wcss[k] <- kmeans.results$tot.withinss
  }
  
  Elbowplot <- plot(1:max_k,wcss,type="b",xlab="Number of clusters", ylab="Within-Cluster Sum of Squares (WCSS)",cex.lab=1.2)
  
}

# Plot the Elbow plot of k-means 
kmeans_Elbow(t(comp_normalized),50,20)

## Determine the number of clusters (Silhouette Method) 
kmeans_silhouette <- function(k){
  km <- kmeans(t(comp_normalized), centers=k, nstart=25)
  ss <- cluster::silhouette(km$cluster, dist(t(comp_normalized)))
  mean(ss[, 3])
}

k <- 2:30
avg_sil <- sapply(k, kmeans_silhouette)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE,cex.lab=1.2)

# another way to visualize silhouette (using factoextra)
factoextra::fviz_nbclust(t(comp_normalized),kmeans,method="silhouette")

## Determine the number of clusters (Gap Statistic Method) 
# compute gap statistic
#' B integer, number of Monte Carlo (bootstrap) samples.
set.seed(123)
# change thr B for parameter optimization (Start with a reasonable range of potential values for B, such as 10, 50, 100, 500, or 1000. )
gap_stat <- cluster::clusGap(t(comp_normalized), FUN=kmeans, nstart=25,
                   K.max=30, B=5)

factoextra::fviz_gap_stat(gap_stat)

## Function to perform k-means clustering 
kmeans_clustering <- function(data,number_of_clusters){
  
  distance_matrix <- as.matrix(amap::Dist(t(data),method="euclidean"))
  # cluster by row, so transpose the matrix to assign cluster to genes
  result <- kmeans(distance_matrix,number_of_clusters,nstart=25,iter.max=1000)
  cluster <- result$cluster
  
  return(cluster)
  
}

## Apply k-means with k =12 (change `number_of_clusters` to 9 and 15 for others result in report)
# perform on expression matrix after dimensionality reduction (comp_normalized)
kmeans_pca_normalized <- kmeans_clustering(comp_normalized,12)
# saveRDS(kmeans_pca_normalized,"kmeans_12_vector.rds")
