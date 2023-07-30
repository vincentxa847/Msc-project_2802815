#### Set up the plasmodium org.db #### [After adjustment]

library("AnnotationForge")
library("tidyr")
library("dplyr")

# Set the paths to the reference files
# Organism
organism <- "P.falciparum_3D7" #change to your organism

# Database
database <- "PlasmoDB" # change to PlasmoDB

# Database release
release <- "release-64" # change to our current release. This is in the header on our sites (should be 64)

# Path to the reference genome
genome <- paste("./rrvis",organism,database,release,"indexes/hisat2/PlasmoDB-64_Pfalciparum3D7_Genome.fasta", sep = "/") #download the reference from Downloads on PlasmoDB and use that path

# Path to the gtf file
gtf <- paste("./rrvis",organism,database,release,"annotation/PlasmoDB-64_Pfalciparum3D7.gff", sep = "/") # download the GFF from Downloads on PlasmoDB and use that path

# Path to GO annotations
go <-  paste("./rrvis",organism,database,release,"annotation/PlasmoDB-CURRENT_Pfalciparum3D7_GO.gaf", sep = "/") # download the GAF from Downloads on PlasmoDB and use that path


# Create the P.falciparum annotation database
# Read in the gene annotations
gff <- as.data.frame(rtracklayer::import(gtf))
gff <- gff[gff$type=="protein_coding_gene",] # 49712 to 5318 
# Replace NA in name column with the gene ID
gff <- gff %>% mutate(Name = ifelse(Name %in% NA, ID, Name))

# Create gene information dataframe
sym <- gff[,c("ID","Name","description")]
colnames(sym) <- c("GID","SYMBOL","GENENAME")
sym <- sym[!duplicated(sym), ]

# Create gene id to chromosome dataframe
chr <- gff[,c("ID","seqnames")]
colnames(chr) <- c("GID","CHROMOSOME")
chr <- chr[!duplicated(chr), ]

# Create GO dataframe
go.data <- read.table(go, header = FALSE, sep = "\t", skip = 1,fill = TRUE)
GO <- go.data[c("V2","V5","V7")]
GO <- GO[GO$V5 != "",]
colnames(GO) <- c("GID","GO","EVIDENCE")
GO <- GO[!duplicated(GO), ]
GO <- subset(GO,grepl("^GO:", GO))

# Get the taxonomy ID from the GO file
taxonomy.id <- go.data$V13[1]
taxonomy.id <- gsub("taxon:","",taxonomy.id)

organism <- "Plasmodium falciparum" #change
genus <- unlist(lapply(strsplit(organism," "),`[[`, 1))
species <- unlist(lapply(strsplit(organism," "),`[[`, 2))

# Create the annotation package and save to current working directory
package <- makeOrgPackage(gene_info = sym,
                          chromosome = chr,
                          go = GO,
                          version = "0.1",
                          maintainer = "Ming-Hsuan Lin <2802815l@student.glasgow.ac.uk>", #change
                          author = "Ming-Hsuan Lin <2802815l@student.glasgow.ac.uk>", #change
                          outputDir = ".",
                          tax_id = taxonomy.id,
                          genus = genus,
                          species = species,
                          goTable = "go")

# Install annoation database
install.packages("./org.Pfalciparum.eg.db", repos=NULL, type = "source") #change

# Load the library
library("org.Pfalciparum.eg.db") #change

# You can now use this like the rrvgo vignette, or in other packages that require an org db.
