library(shiny)
library(ggvis)
library(reshape2) # for GOenrichment function 
library(topGO) # for GOenrichment function
library(ComplexHeatmap) # for interactive heatmap
library(InteractiveComplexHeatmap) # for interactive heatmap
# shiny [version 1.7.4]
# ggvis [version 0.4.8]
# InteractiveComplexHeatmap [version 1.2.0]
# dplyr [version 1.1.0]
# shinycssloaders [version 1.0.0.9011]
# DT [version 0.28]


#### Set up data and function ####

## GO enrichment 
geneId2GO <- readRDS(file="geneId2GO.rds")
GOenrichment <- readRDS(file="GOenrichment.rds")

data <- readRDS(file="comp_normalized.rds") # first 200 PCs of normalized data 
original_matrix <- readRDS("Matrix_normalized_excluded_removeanti.rds")

## Gene information table
genes_description <- read.csv("5576_genes_with_description.txt",sep = "\t")
genes_description <- genes_description[!duplicated(genes_description[,1]),] # remove row with duplicated GENE.id

## UMAP
umap_coordinate <- readRDS(file="umap_coordinate.rds")
umap_coordinate <- as.data.frame(umap_coordinate$layout)
umap_coordinate$gene <- row.names(umap_coordinate)

## knn graph
knn_graph <- readRDS(file="knn_graph_louvain_jaccard_threshold.rds")



#### Define UI ####

ui <- fluidPage(
  
  titlePanel("Louvain Clustering"), # Add title 
  
  fluidRow(
    column(3,
           wellPanel(
             selectInput("method",label="Clustering method :",
                         choices="louvain clustering"),
             numericInput("resolution", label="Resolution :", 
                          value=0.6,min=0.1,max=10)
           ),
           wellPanel(
             numericInput("Filter", label="Select the cluster. Essential for GO enrichment (0 is select all) :",
                          value = 0, min=1, max=20),
             textInput("gene_filter", "Gene(s) to search", "Using `,` to seperate genes"),
             sliderInput("X_axis", "Range of UMAP 1 to apply",
                         -7.5,7.5,c(-7.5,7.5), step = 0.1),
             sliderInput("Y_axis", "Range of UMAP 2 to apply", -7.5,7.5,c(-7.5,7.5), step = 0.1),
             actionButton("submitbutton", "Submit", icon("paper-plane"))
           ),
    ),
    column(9,
           mainPanel(
             tabsetPanel(
               tabPanel("UMAP Plot", 
                        ggvisOutput("UMAPplot"), 
                        DT::DTOutput("Genetable")),
               
               tabPanel("Heatmap", 
                        shinycssloaders::withSpinner(InteractiveComplexHeatmapOutput()),
                        # Create interactive complex heatmap take a lot of time, so use another active button
                        actionButton("heatmapbutton", "Create heatmap", class="btn btn-outline-secondary")),
               
               tabPanel("GO enrichment", 
                        br(),
                        br(),
                        DT::DTOutput("GOtable"),
                        actionButton("GObutton", "Run GO enrichment", class="btn btn-outline-secondary")),
               
               tabPanel("GSEA", 
                        br(),
                        selectInput("Sample",label="Choose a sample for GSEA:",
                                    # set the width = "100%" to make sure input display the whole datasets name
                                    choices=row.names(original_matrix),width="100%"), 
                        actionButton("GSEAsubmit", "Submit", icon("paper-plane")),
                        br(),
                        br(),
                        DT::DTOutput("GSEAtable"))
             )
           )
    )
  ) 
)




#### Define server ####

server <- function(input,output,session) {
  
  
  ## Set up the information that user apply
  clusters <- reactive({
    
    # Louvain algorithm 
    # knn_igraph here using k = 15 after k optimization 
    louvain_resolution <- input$resolution
    method <- input$method
    
    if (method == "louvain clustering"){
      clusterlouvain <- igraph::cluster_louvain(knn_graph ,resolution = louvain_resolution ) 
      louvain_cluster <- clusterlouvain$membership
      names(louvain_cluster) <- colnames(original_matrix)
    }
    cluster_to_use <- louvain_cluster
  })
  
  ## Function for generating tooltip text (being used in interactive UMAP plot -- line ??)
  gene_tooltip <- function(x) { # the input (x) from ggvis is V1,V2 and gene, can use $ to specify
    if (is.null(x)) return(NULL) # without get the input from ggvis
    if (is.null(x$gene)) return(NULL)
    
    cluster <- isolate(clusters()) # clusters should be named vector
    
    # some genes have more than one transcript, so using gene id rather row name here
    gene_description <- genes_description[genes_description$Gene.ID == as.character(x$gene),3]
    result <- cluster[as.character(x$gene)] 
    
    # Show the gene name and cluster belong to (cluster information according to user select)
    paste0("Gene name :", "<b>",x$gene,"</b><br>",
           "Description :","<b>",gene_description,"</b><br>",
           "Cluster :","<b>",result, "</b><br>")
    
  }
  
  ## A reactive subset of umap_coordinate *was designed to zooming the plot*
  # filter function from dplyer
  selected <- reactive({ 
    umap_coordinate_selected <- umap_coordinate %>% 
      dplyr::filter(V1 > input$X_axis[1] & V1 < input$X_axis[2]) %>%
      dplyr::filter(V2 > input$Y_axis[1] & V2 < input$Y_axis[2]) 
    
    umap_coordinate_selected
  })
  
  ## Set up GO enrichment
  observeEvent(input$GObutton,{
    
    req(input$Filter) # exit the observeEvent if user not input any nmuber
    
    # GO enrichment function take two arguments, data (named int) and which_cluster (numeric)
    data <- clusters()
    which_cluster <- input$Filter
    
    shinycssloaders::showPageSpinner() # loading animation
    if (which_cluster == 0) {
      # pass, when user not specify the cluster
      GOresult <- data.frame()
    } else {
      GOresult <- GOenrichment(data,which_cluster) # return a table
    }
    shinycssloaders::hidePageSpinner() # loading animation
    
    output$GOtable <- DT::renderDT(server=FALSE,{
      
      DT::datatable(
        data=GOresult,
        extensions='Buttons',
        options=list(scrollY = "300px",
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'))
      )
    })
  })
  
  ## Set up GSEA table
  observeEvent(req(input$Sample, input$method, input$GSEAsubmit),{ # observeEvent with multiple conditions
    
    clusters <- isolate(clusters()) # clusters should be named vector
    sample_chosed <- input$Sample
    
    # pathways -- list of the genes in each cluster
    # create a vector with empty list
    gseaList <- vector("list", length(unique(as.character(clusters))))
    
    for (cluster in unique(as.character(clusters))) {
      geneList <- names(clusters[which(clusters == cluster)])
      # push the geneList into list (in order)
      gseaList[[paste0('cluster', cluster)]] <- as.character(geneList)
    }
    
    # Run GSEA
    shinycssloaders::showPageSpinner()  # loading animation
    
    ranks <- original_matrix[sample_chosed,]
    print(sample_chosed)
    # names(ranks) <- df$GeneId
    ranks <- ranks[order(ranks)]
    fgseaRes <- fgsea::fgsea(gseaList, ranks)
    # order by p value
    fgseaRes <- fgseaRes[order(pval)] 
    
    shinycssloaders::hidePageSpinner()  # loading animation

    GSEA_table <- fgseaRes[,c(1:7)] # exclude leading edge
    
    output$GSEAtable <- DT::renderDT(server=FALSE,{
      
      DT::datatable(
        data=GSEA_table,
        extensions='Buttons',
        options=list(scrollY = "300px",
                       scrollX= TRUE,
                       autoWidth = TRUE,
                       order = list(list(6, 'desc')), # use NES to rank
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'))
      )
    })
  })
  
  
  
  ## Create interactive UMAP plot and corresponding gene table
  observeEvent(input$submitbutton,{
    
    req(input$Filter) # exit the observeEvent if user not input any nmuber 
    
    # using *observeEvent* here so not need to dynamically change the output and use *isolate* 
    coordinate <- isolate(selected()) 
    cluster <- isolate(clusters()) 
    gene_table <- genes_description
    
    # Set up the filter of cluster
    cluster_choiced <- input$Filter
    # 0 is select all
    if (cluster_choiced != 0) {
      clusterAFTERfilter <-  cluster[cluster == cluster_choiced]
      
      #' in a situation that coordinate have been updated, but the step to update cluster based on
      #' coordinate is later on line 162 *cluster <- cluster[names(cluster) %in% row.names(coordinate)]* 
      #' and the  *clusterAFTERfilter* to update coordinate is based on cluster, so this is irrational
      #' therefore update the *clusterAFTERfilter* here
      clusterAFTERfilter <- clusterAFTERfilter[names(clusterAFTERfilter) %in% row.names(coordinate)]
      
      coordinate <- coordinate[names(clusterAFTERfilter),] 
      
      gene_table <- subset(gene_table,gene_table[,"Gene.ID"] %in% names(clusterAFTERfilter))
      gene_table$cluster <- clusterAFTERfilter
    } else {
      gene_table$cluster <- cluster
    }
    
    # Set up the filter of gene(s)
    gene_selected <- input$gene_filter
    if (gene_selected != "Using `,` to seperate genes" & gene_selected != "") {
      # multiple genes input
      if (grepl(",",gene_selected)){
        gene_selected <- as.character(unlist(strsplit(gene_selected,",")))
      }
      
      
      cluster <- cluster[gene_selected]
      coordinate <- coordinate[gene_selected,]
      gene_table <- subset(gene_table,gene_table[,"Gene.ID"] %in% gene_selected)
      gene_table$cluster <- cluster
    }
    
    
    # exclude the Input.ID
    gene_table <- gene_table[,-5]
    
    output$Genetable <- DT::renderDT(server=FALSE,{
      
      DT::datatable(
        data=gene_table, 
        extensions='Buttons',
        options=list(scrollY = "300px",
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'))
      )
    })
    
    # when updating the coordinate (zooming), update the cluster too
    cluster <- cluster[names(cluster) %in% row.names(coordinate)]
    
    
    main_plot <- coordinate %>%
      ggvis(x= ~V1,y= ~V2, key := ~gene) %>%
      layer_points(fill = ~factor(cluster)) %>% 
      add_tooltip(gene_tooltip,"hover") %>%
      add_axis("x", title = "UMAP 1") %>%
      add_axis("y", title = "UMAP 2")  %>%
      add_legend("fill", title = "Cluster") %>%
      scale_nominal("fill", range = c("#8f5362","#BF9189","#9C7754","#e0a981","#ecd09c","#d4daa1","#a3c8a4","#79b4a0",
                                      "#6888a5","#369F2D","#FAC7B3","#EE4431",
                                      "#1663A9","#B4B4D5","#8481BA","#B0A875"))%>% 
      set_options(height = 700, width = 900)
    
    ## Render the plot
    main_plot %>% bind_shiny("UMAPplot")
    
  })
  
  ## Create interactive heatmap
  observeEvent(input$heatmapbutton,{
    
    column_split <- clusters() # named int
    # cluster_choiced <- input$Filter
    
    # Define the color gradient
    my_palette <- colorRampPalette(c("#C25160", "#5AA4AE","#806D9A","#4F794A","#06436F"))
    
    # plot on matrix before PCA
    dataTOplot <- original_matrix
    
    shinycssloaders::showPageSpinner()  # loading animation
    ht1 <- ComplexHeatmap::Heatmap(dataTOplot, name = "TPM",
                                  column_split = column_split,show_column_dend = FALSE,
                                  cluster_rows = TRUE, show_row_dend = FALSE,
                                  show_row_names = FALSE,row_names_gp = gpar(fontsize = 10),
                                  show_column_names= FALSE,
                                  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = my_palette(20)))))
    ht1 <- draw(ht1)
    
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input,output,session,ht_list = ht1)
    shinycssloaders::hidePageSpinner()  # loading animation
  })
  
  
  
}
shinyApp(ui, server)

