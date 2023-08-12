## change data structure
clinical_isolate <- readRDS("./Datasets_directory/clinical_isolate.rds")
oocyst <- readRDS("./Datasets_directory/oocyst.rds")
gametocyte <- readRDS("./Datasets_directory/gametocyte.rds")
schizont <- readRDS("./Datasets_directory/schizont.rds")
troph <- readRDS("./Datasets_directory/troph.rds")
ring <- readRDS("./Datasets_directory/ring.rds")
Asexual_blood_stage <- readRDS("./Datasets_directory/Asexual_blood_stage.rds")
merozoite <- readRDS("./Datasets_directory/merozoite.rds")
sporozoite <- readRDS("./Datasets_directory/sporozoite.rds")

# Set the source directory path to where the GSEA files located
source_directory <- "./GSEA_hclust_12"


# List all the files in the source directory
files <- list.files(path = source_directory, full.names = TRUE)

# Loop through each file
for (file in files) {
  # Extract the file name without the path
  file_name <- basename(file)
  
  # Check if the file belongs to the "ring" category
  if ( file_name %in% ring) {
    # Set the destination subdirectory for "ring" files
    destination_subdirectory <- paste0(source_directory,"/ring")
  }
  # Check if the file belongs to the "troph" category
  else if ( file_name %in% troph) {
    # Set the destination subdirectory for "troph" files
    destination_subdirectory <- paste0(source_directory,"/troph")
  }
  # Check if the file belongs to the "schizont" category
  else if (file_name %in% schizont) {
    # Set the destination subdirectory for "schizont" files
    destination_subdirectory <- paste0(source_directory,"/schizont")
  }
  # Check if the file belongs to the "gametocyte" category
  else if (file_name %in% gametocyte) {
    # Set the destination subdirectory for "gametocyte" files
    destination_subdirectory <- paste0(source_directory,"/gametocyte")
  }
  
  # Check if the file belongs to the "asexual blood stage" category
  else if (file_name %in% Asexual_blood_stage) {
    # Set the destination subdirectory for "merozoite" files
    destination_subdirectory <- paste0(source_directory,"/Asexual_blood_stage")
  }
  # Check if the file belongs to the "clinical isolate" category
  else if (file_name %in% clinical_isolate) {
    # Set the destination subdirectory for "schizont" files
    destination_subdirectory <- paste0(source_directory,"/clinical_isolate")
  }
  # Check if the file belongs to the "oocyst" category
  else if (file_name %in% oocyst) {
    # Set the destination subdirectory for "oocyst" files
    destination_subdirectory <- paste0(source_directory,"/oocyst")
  }
  # Check if the file belongs to the "merozoite" category
  else if (file_name %in% merozoite) {
    # Set the destination subdirectory for "merozoite" files
    destination_subdirectory <- paste0(source_directory,"/merozoite")
  }
  # Check if the file belongs to the "sporozoite" category
  else if (file_name %in% sporozoite) {
    # Set the destination subdirectory for "sporozoite" files
    destination_subdirectory <- paste0(source_directory,"/sporozoite")
  }
  else {
    destination_subdirectory <- paste0(source_directory)
  }
  
  # Create the destination subdirectory if it doesn't exist
  dir.create(destination_subdirectory, recursive = TRUE, showWarnings = FALSE)
  
  # Construct the new file path with the destination subdirectory
  new_file_path <- file.path(destination_subdirectory, file_name)
  
  # Move the file to the destination subdirectory
  file.rename(from = file, to = new_file_path)
}


