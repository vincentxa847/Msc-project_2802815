#### DBSCAN ####
# nndist from spatstat.geom [version 3.1-0]
# dbscan [version 1.1-11]

comp_normalized <- readRDS(file="comp_normalized.rds")

## Parameters optimization
# choose minPts 400 according to rule of thumb
# using likes elbow method to choose optimal eps


# `nndist` Computes the distance from each point to its nearest neighbour in a point pattern.
kdist_DBSCAN <- spatstat.geom::nndist(t(comp_normalized), k=400, by=NULL, method="C", metric=NULL)
sorted_kdist_DBSCAN <- sort(kdist_DBSCAN, decreasing=FALSE)

sorted_kdist_DBSCAN_matrix <- data.frame(k=seq(from=1, to=length(sorted_kdist_DBSCAN)), distance=sorted_kdist_DBSCAN)
ggplot(sorted_kdist_DBSCAN_matrix) +
  geom_point(aes(x = k, y = distance)) +
  labs(x = "datapoint", y = "distance") + theme_bw()

## Perform DBSCAN
# from parameter optimization, choosing `eps = 5 , minPts = 400` here
# genes and their corresponding clusters are stored in dbscan_result$cluster
dbscan_result <- dbscan::dbscan(t(comp_normalized),eps=5,minPts=400)

# change the standard of dense by changing minPts to 3
dbscan_result <- dbscan::dbscan(t(comp_normalized),eps=5,minPts=3)
