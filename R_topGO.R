source("https://bioconductor.org/biocLite.R")

library("topGO")
library("biomaRt")

args = commandArgs(trailingOnly=TRUE)


mapping_db <- args[5]
path_results <- args[4]
targets_file <- args[2]
universe_file <- args[3]
ID <- args[1]

universe <- read.table(universe_file, header = FALSE, sep = "\t")
targets <- read.table(targets_file, header = FALSE, sep = "\t")

allgeneNames <- universe$V1
geneNames<-allgeneNames[universe$V1]

geneList <- factor(as.integer(geneNames %in% targets$V1))
names(geneList) <- geneNames

# # Run topGO with Cellular Component as the ontology

if (mapping_db == "org.At.tair.db") {
	CC_GOdata <- new("topGOdata",ontology = "CC", allGenes = geneList, geneSelectionFun = function(x)x, annot = annFUN.org, mapping = mapping_db)
} else {
	CC_GOdata <- new("topGOdata",ontology = "CC", allGenes = geneList, annot=annFUN.org, mapping=mapping_db, ID="Ensembl")
}


CC_resultFisher <- runTest(CC_GOdata , algorithm = "classic", statistic = "fisher")
CC_resultKS <- runTest(CC_GOdata , algorithm = "classic", statistic = "ks")
CC_resultKS.elim <- runTest(CC_GOdata , algorithm = "elim", statistic = "ks")

CC_allRes <- GenTable(CC_GOdata, classicFisher = CC_resultFisher, classicKS = CC_resultKS, elimKS = CC_resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)


# Run topGO with Biological Process as the ontology
if (mapping_db == "org.At.tair.db") {
	BP_GOdata <- new("topGOdata",ontology = "BP", allGenes = geneList, geneSelectionFun = function(x)x, annot=annFUN.org, mapping=mapping_db)
} else {
	BP_GOdata <- new("topGOdata",ontology = "BP", allGenes = geneList, annot=annFUN.org, mapping=mapping_db, ID="Ensembl")
}

BP_resultFisher <- runTest(BP_GOdata , algorithm = "classic", statistic = "fisher")
BP_resultKS <- runTest(BP_GOdata , algorithm = "classic", statistic = "ks")
BP_resultKS.elim <- runTest(BP_GOdata , algorithm = "elim", statistic = "ks")

BP_allRes <- GenTable(BP_GOdata, classicFisher = BP_resultFisher, classicKS = BP_resultKS, elimKS = BP_resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

# Write the tables to files

# Define the name of the table
name_file_topGO_CC <- paste('TableTopGO_CC_',ID,'.txt', sep='')
name_file_topGO_BP <- paste('TableTopGO_BP_',ID,'.txt', sep='')

# Define the file paths
path_to_table_topGO_CC <- paste(path_results, name_file_topGO_CC, sep='/')
path_to_table_topGO_BP <- paste(path_results, name_file_topGO_BP, sep='/')

# Write the output topGO tables
write.table(CC_allRes, path_to_table_topGO_CC, sep="\t", quote=FALSE)
write.table(BP_allRes, path_to_table_topGO_BP, sep="\t", quote=FALSE)











