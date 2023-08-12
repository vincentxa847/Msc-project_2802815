## From the view of data sets
# Fetch all available sets
clinical_isolate <- readRDS("./Datasets_directory/clinical_isolate.rds")
oocyst <- readRDS("./Datasets_directory/oocyst.rds")
gametocyte <- readRDS("./Datasets_directory/gametocyte.rds")
schizont <- readRDS("./Datasets_directory/schizont.rds")
troph <- readRDS("./Datasets_directory/troph.rds")
ring <- readRDS("./Datasets_directory/ring.rds")
Asexual_blood_stage <- readRDS("./Datasets_directory/Asexual_blood_stage.rds")
merozoite <- readRDS("./Datasets_directory/merozoite.rds")
sporozoite <- readRDS("./Datasets_directory/sporozoite.rds")



list.dirs()
path = "./GSEA_hclust_12/troph" # change here to the GSEA directory and interesting stage
setwd(path)
getwd()

test = data.frame()
# "troph" here represents files corresponding to this stage
for (file in unique(as.character(troph))) { # change here for other interesting stage
  table <- openxlsx::read.xlsx(file)
  
  # set the p value threshold and pick the positive enrichment
  table <- subset(table, padj < 0.05 & NES >0)
  
  # row.names(table) = as.character(file)
  test = rbind(test,table)
  
}
test$pathway = factor(test$pathway)
library(ggplot2)
ggplot(test,aes(x=factor(pathway),y=NES, colour=padj))+ geom_point() + theme_minimal() +
	# change the font size for x axis label and title
	theme(axis.text.x = element_text(size = 14),plot.title = element_text(size = 18)) +
    stat_summary(
    fun.data = function(x) data.frame(y = max(x), label = length(x)),
    geom = "text",
    position = position_nudge(y = 0.1),
    size = 5,
    show.legend = FALSE
    ) + labs(x="",y="NES",title = "Schizont  (66 samples)") + # change here for visualizing other stages 
  scale_colour_continuous(low="red", high="blue",guide=guide_colorbar(reverse=TRUE))

## From the view of cluster [prefer this one]
# change here to the GSEA directory 
path = "C:/Users/user/Desktop/2802815_Msc_Project/GSEA/GSEA_hclust_12"
setwd(path)
all_files = list.files(path=path,recursive = TRUE) 

GSEA_visualization <- function(cluster, showDataset=30, font.size=4) {
  
  cluster_ <- data.frame()
  cluster <- paste0("cluster",as.character(cluster))
  
  for (file in unique(as.character(all_files))) {
    table <- openxlsx::read.xlsx(file)
    
    # set the p value threshold and pick the positive enrichment
    table <- subset(table, padj <0.001 & NES >0)
    if (cluster %in% table$pathway) {
      # create a table stroing datasets with cluster 1 significant enriched 
      print(file)
      test = subset(table,table$pathway == cluster)
      print(test[,c(1:7)])
      test$dataset = as.character(gsub("_GSEA.xlsx","",file))
      cluster_ = rbind(cluster_,test)
    }
  }
  
  cluster_idx = order(cluster_[["NES"]],decreasing = TRUE)
  cluster_idx = cluster_[cluster_idx,]
  cluster_ = cluster_idx[c(1:showDataset),]
  
  
  # get the NES order 
  idx = order(cluster_[["NES"]],decreasing = TRUE)
  # change the dataset odrer so the dataset in y axis can order in decreasing order [change factor level, which is invisual in matrix]
  cluster_$dataset <- factor(cluster_$dataset,
                              levels=rev(unique(cluster_$dataset[idx]))) 
  
  ggplot(cluster_, aes(x=NES, y=dataset, colour=padj)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue",
                           guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(paste0(cluster," (",as.character(length(idx))," samples",")")) +
    theme_minimal()+
    theme(axis.text.y = element_text(size = font.size)) + # theme work only if behind theme_minimal
    guides(size  = guide_legend(order = 1),
           color = guide_colorbar(order = 3))
}

GSEA_visualization(12, showDataset=30, font.size = 10) # change here for different cluster visualization 
