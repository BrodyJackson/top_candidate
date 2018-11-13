

install("topCandidate")
library(topCandidate)

dataList <- ortholog_import("all_orthologs_pine_spruce.txt","one_to_one_orthologs.txt")
allOrthologs <- data.frame(dataList[[1]])
oneToOne <- data.frame(dataList[[2]])

pine <- species_data_import("all_pine_data_pvalues.txt", TRUE)
overall_pine <- data.frame(dataList[[1]])
pine_tcontig <- data.frame(dataList[[2]])

spruce <- species_data_import("all_spruce_data_pvalues.txt", TRUE)
overall_spruce <- data.frame(dataList[[1]])
spruce_tcontig <- data.frame(dataList[[2]])

ortholog_cleaning(pine_tcontig, spruce_tcontig, 3)



