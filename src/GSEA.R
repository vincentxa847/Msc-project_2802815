## Create sub directory to save the GSEA result table 
dir.create("./GSEA_new_kmeans_11")
setwd("./GSEA_new_kmeans_11")
getwd() 

## cluster
clusters <- kmeans_pca_normalized_11
clusters <- as.character(clusters)
names(clusters) <- names(kmeans_pca_normalized_11)
clusters

## expression data
expression <- Matrix_normalized_excluded_removeanti

## pathways -- list of the genes in each cluster
# create a vector with empty list
gseaList <- vector("list", length(unique(as.character(clusters))))
for (cluster in unique(as.character(clusters))) {
  geneList <- names(clusters[which(clusters == cluster)])
  # push the geneList into list (in order)
  gseaList[[paste0('cluster', cluster)]] <- as.character(geneList)
}
gseaList


## iterate the fgsea for each sample
for (sample in row.names(Matrix_normalized_excluded_removeanti)) {
  print(sample)
  ranks <- Matrix_normalized_excluded_removeanti[sample,]
  # names(ranks) <- df$GeneId
  ranks <- ranks[order(ranks)]
  fgseaRes <- fgsea::fgsea(gseaList, ranks)
  # print(fgseaRes)
  # order by p value
  fgseaRes <- fgseaRes[order(pval)] 
  file_name = paste0(sample,"_GSEA",".xlsx")
    openxlsx::write.xlsx(fgseaRes,file = file_name,rownames = TRUE)
}
  
setwd("..") # leave subdirectory
getwd()
