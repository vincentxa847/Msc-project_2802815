#### Loading TPM file ####

# load the srr,err and their corresponding condition 
SRR_and_ERR_for_index <- readRDS(file="SRRandERRTable.rds") # from parsing xml file using Python and R
condition <- SRR_and_ERR_for_index[,1]
srr <- SRR_and_ERR_for_index[,2]

## function to access file path
file_path <- function(resource) {
  directory <- paste("./Project_Data/",resource,sep="")
  sample <- list.files(directory)
  path <- paste(directory,sample,sep="/")
  return(path)
}

## Dataset of VEuPathPipline and ebiPipline TPM
TPM_VEuPathPipline <- file_path("VeuPathPipeline_TPM_output")
TPM_ebiPipline <- file_path("ebiPipline_TPM_output")

## function to create matrix
create_matrix <- function(files) {
  data <- c()
  column_names <- c()
  for (file in files) {
    # get the experiment group and sample information (as column of matrix) and gene name (as row of matrix) 
    string <- strsplit(file, "/")
    
    experimentANDsample <- string[[1]][4]
    # run the check to replace srr or err number to condition
    for (check in seq_along(srr)) {
      k <- srr[check]
      z <- condition[check]
      if (grepl(k,experimentANDsample)){
        experimentANDsample <- gsub(k,z,experimentANDsample)
      } 
    } 
    experimentANDsample <- gsub("genes.htseq-union.","",experimentANDsample)
    
    column_name <- experimentANDsample 
    column_name <- gsub("TPM.","",column_name)
    column_name <- gsub(".TPM","",column_name)
    column_name <- gsub(".counts","",column_name)
    
    column_names <- c(column_names,column_name)
    
    # open file
    file <- read.csv(file,header = FALSE,row.names = 1,col.names = c("1",experimentANDsample),skip = 1,sep = "\t")
    gene.name <- row.names(file)
    
    # concatenate
    data <- c(data,file)
  }
  
  # Create matrix
  Matrix <- matrix(unlist(data),nrow = 5720)
  rownames(Matrix) <- gene.name
  colnames(Matrix) <- column_names
  
  # Transpose the matrix to make genes as columns (genes here is like cells in expression matrix)
  Matrix <- t(Matrix)
  
  return(Matrix)
}

#### Combine both ebi and VEuPathDB ####

Matrix_VEuPath <- create_matrix(TPM_VEuPathPipline) # 36
Matrix_ebi <- create_matrix(TPM_ebiPipline) # 570
Matrix <- rbind(Matrix_VEuPath,Matrix_ebi) # 606

## Perform z score normalization
Matrix_normalized <- scale(Matrix,center = TRUE,scale = TRUE)

# which(colSums(is.nan(Matrix_normalized))>0)
# the combined matrix have only one gene with NaN *PF3D7_1200610*

# Exclude genes with NaN 
Matrix_normalized <- Matrix_normalized[,-c(3532)]

#### Remove rRNA and tRNA genes from matrix ####

rRNA_tRNA <- read.csv("rRNAandtRNA_3D7.csv")[,1] # length 143

#' names function can only be applied to data.frame. to get the column names
#' %in% operator create TRUE/FALSE vector to indicate which column to remove
#' if the list contain name of row or column, then -c() cannot be used, rather need to use %in%
Matrix_normalized_excluded <- Matrix_normalized[, !(names(as.data.frame(Matrix_normalized)) %in% rRNA_tRNA)]


#### Function to remove antisense strand ####
removeantistrand <- function(matrix_to_used) {
  
  matrix <- matrix_to_used
  ## Exclude antisense strsnd
  # `grepl` is to test if rownames(Matrix) have firststrand and return logical vector
  firststrands <- row.names(matrix[grepl("firststrand",rownames(matrix)),])

  strandToRemove <- c()
  for (i in firststrands) {
    # first strand name and second strand name 
    firststrand <- i
    secondstrand <- gsub("firststrand","secondstrand",i)
    
    # sum the TPM value of each strand across genes
    first <- sum(matrix[firststrand,])
    second <- sum(matrix[secondstrand,])
    
    if (first > second) {
      strandToRemove <- c(strandToRemove,secondstrand)
    }
    
    else if (second > first){
      strandToRemove <- c(strandToRemove,firststrand)
    }
    
  }
    
    Matrix_after_processed <- matrix[ !(rownames(as.data.frame(matrix)) %in% strandToRemove),]
    
    return(Matrix_after_processed)
    
}

Matrix_normalized_excluded_removeanti <-  removeantistrand(Matrix_normalized_excluded)
# saveRDS(Matrix_normalized_excluded_removeanti,file="Matrix_normalized_excluded_removeanti.rds")
