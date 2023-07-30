#### Perform PCA on TPM matrix ####
# prcomp from base package stats [version 4.1.3]

#' prcomp expects to see components (samples) by row and dimensions (genes) by column,
#' turn off scale here

Matrix_normalized_excluded_removeanti <- readRDS(file="Matrix_normalized_excluded_removeanti.rds")
pca_normalized <- prcomp(t(Matrix_normalized_excluded_removeanti))

PCA_ElbowPlot <- function(pca) {
  # Calculate the variance that each component explained
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  # Cumulative proportion of variance explained
  cumulative_variance <- cumsum(variance_explained)
  # Elbow data (principal components and corresponding cumulative proportion of variance explained)
  elbow.data <- data.frame(principal_components = 1:length(cumulative_variance),
                          cumulative_variance = cumulative_variance)
  
  # plot
  plot(elbow.data$principal_components, elbow.data$cumulative_variance,
       type = "l", xlab = "Number of principal components",
       ylab = "Cumulative variance explained")
  points(elbow.data$principal_components, elbow.data$cumulative_variance, type = "p")
  
  # Add vertical dotted line at x = 200
  abline(v = 200, lty = "dotted")
}

elbowplot_normalized <- PCA_ElbowPlot(pca_normalized)

## Take the first 200 PCs of normalized data 
comp_normalized <- data.frame(t(pca_normalized$x[,1:200]))
# saveRDS(comp_normalized, file = "comp_normalized.rds")
