library("reshape2")
library("topGO")

geneId2GO <- file("TopGO/GenesByTaxon_Summary.txt") 
geneId2GO <- readLines(geneId2GO)
geneId2GO <- gsub(';', '\t', geneId2GO)
geneId2GO <- strsplit(geneId2GO, '\t')
n <- sapply(geneId2GO, length)
seq.max <- seq_len(max(n))
geneId2GO <- data.frame(t(sapply(geneId2GO, "[", i = seq.max)))

geneId2GO[geneId2GO=="N/A"] <- NA
geneId2GO <- geneId2GO[,c(1,3:ncol(geneId2GO))]

geneId2GO <- geneId2GO[-1,]

geneId2GO <-melt(geneId2GO, id.vars='X1', na.rm=T)

geneId2GO <- geneId2GO[,c(1,3)]

colnames(geneId2GO) <- c("geneId", "value")
geneId2GO[sapply(geneId2GO, is.factor)] <- lapply(geneId2GO[sapply(geneId2GO, is.factor)],
                                                  as.character)
geneId2GO <-melt(geneId2GO, id.vars='geneId', na.rm=T)
geneId2GO <- geneId2GO[,c(1,3)]


geneId2GO <- by(geneId2GO$value,
                geneId2GO$geneId,
                function(x) as.character(x))
head(geneId2GO)

# saveRDS(geneId2GO, file = "geneId2GO.rds")

GOenrichment = function(data,which_cluster) {
  ## Defining your list of genes of interest, and the 'gene universe' to compare it to
  geneUniverse <- names(data)
  genesOfInterest <- names(data[data == which_cluster])
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # genelist must be 2 level, level 0 background and level 1 gene of interest
  head(geneList)  

  ## create a topGO object
  # geneList is a list of genes. In your case, this would be a list of the genes in a cluster
  # ontology can be BP (Biological Process), MF (Molecular Function) or CC (Cellular Compartment)
  BPenrichment <- new("topGOdata", description=gsub("XX",as.character(which_cluster),"Cluster XX Biological Process Enrichment Analysis"), 
                      ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneId2GO)

  # show genes in your list that have BP terms and can contribute to the analysis
  # sigGenes(BPenrichment)

  ## do the enrichment analysis using Fisher test with weight01 pruning algorithm and show the results
  BPFisherResults <- runTest(BPenrichment, statistic="fisher")
  # BPFisherResults

  ## Make a list of the enriched terms. Here, I have opted to only display GO terms with P < 0.05 in the table
  BPenriched <- GenTable(BPenrichment, Fisher=BPFisherResults, orderBy="Fisher", ranksOf=BPFisherResults, topNodes=20)
  BPsignificant <- BPenriched[ which(as.numeric(BPenriched$Fisher) <= 0.05), ]
  # BPsignificant
  
  # write out the table
  # subdirectory = 'TopGO'
  # file_name = paste(deparse(substitute(data)),as.character(which_cluster),sep = "_")
  # write.table(BPsignificant,file.path(subdirectory,file_name),sep = "\t")
  
	# return a table
  return(BPsignificant)
}

# saveRDS(GOenrichment, file = "GOenrichment.rds")

# Perform GO analysis
clustering_result <- readRDS(file="../Louvain_clustering/louvain_0.6_new.rds")
which_cluster <- 1
GO <- GOenrichment(clustering_result,which_cluster)

