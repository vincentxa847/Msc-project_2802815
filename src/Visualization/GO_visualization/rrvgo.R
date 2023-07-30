# rrvgo [version 1.6.0]

library("org.Pfalciparum.eg.db")
data <- kmeans_pca_normalized_18
which_cluster <- 1
## Defining your list of genes of interest, and the 'gene universe' to compare it to
geneUniverse <- names(data)
genesOfInterest <- names(data[data == which_cluster])
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# genelist must be 2 level, level 0 background and level 1 gene of interest
head(geneList)  
BPenrichment <- new("topGOdata", description= gsub("XX",as.character(which_cluster),"Cluster XX Biological Process Enrichment Analysis"), 
                    ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneId2GO)
# show genes in your list that have BP terms and can contribute to the analysis
sigGenes(BPenrichment)

## do the enrichment analysis using Fisher test with weight01 pruning algorithm and show the results
BPFisherResults <- runTest(BPenrichment, statistic="fisher")
BPFisherResults
BPenriched <- GenTable(BPenrichment, Fisher=BPFisherResults, orderBy="Fisher", ranksOf=BPFisherResults, topNodes=20)
BPsignificant <- BPenriched[ which(as.numeric(BPenriched$Fisher) <= 0.05), ]
# BPsignificant


simMatrix <- rrvgo::calculateSimMatrix(BPsignificant$GO.ID, orgdb="org.Pfalciparum.eg.db", ont="BP", method="Rel")


scores <- setNames(BPsignificant$Fisher, BPsignificant$GO.ID)

reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Pfalciparum.eg.db")

# Heatmap is suit for the report
rrvgo::heatmapPlot(simMatrix, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=6)
