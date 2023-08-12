# get.knn from [FNN version 1.1.3.1]
# igraph version 1.4.2

comp_normalized <- readRDS(file="../PCA/comp_normalized.rds")
Matrix_normalized_excluded_removeanti <- readRDS(file="../Data_Preparation/Matrix_normalized_excluded_removeanti.rds")

x <- t(comp_normalized)
k <-15
n <- nrow(x)

# adapt from cccp package nng function
# using get.knn from FNN 
dx <- FNN::get.knn(x,k=k,algorithm="cover_tree")

edges <- matrix(unlist(sapply(1:nrow(x),function(i) {
  rbind(rep(i,k),dx$nn.index[i,])
})),nrow=2) # length 83640 (two rows, store the 2 nodes of each cluster)

edges_jaccard_0.05 <- edges

for (edge in 1:ncol(edges)) {
  first_node = edges[1,edge]
  second_node = edges[2,edge]
  jaccard_index = length(intersect(dx$nn.index[first_node,],dx$nn.index[second_node,]))/
    length(unique(c(dx$nn.index[first_node,],dx$nn.index[second_node,])))
  if (jaccard_index < 0.05){
    edges_jaccard_0.05 = edges_jaccard_0.05[,-edge]
  }
  
}
length(edges_jaccard_0.05[1,]) # 73334

out <- igraph::make_graph(edges=edges,n=n,directed=FALSE)
saveRDS(out, file="knn_graph_louvain_jaccard_threshold.rds")
# igraph::plot.igraph(out, vertex.label = NA,vertex.size=5)


clusterlouvain <- igraph::cluster_louvain(out,resolution = 0.6) 
louvain_cluster <- clusterlouvain$membership
names(louvain_cluster) <- colnames(Matrix_normalized_excluded_removeanti)
length(unique(clusterlouvain$membership))

# UMAP_visualization(umap_coordinate,louvain_cluster)

saveRDS(louvain_cluster, file="louvain_0.6_new.rds") # 12 clusters
