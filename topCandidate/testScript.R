

install("topCandidate")
library(topCandidate)

dataList <- ortholog_import("all_orthologs_pine_spruce.txt","one_to_one_orthologs.txt")
allOrthologs <- data.frame(dataList[[1]])
oneToOne <- data.frame(dataList[[2]])

pine <- species_data_import("all_pine_data_pvalues.txt", TRUE)
overall_pine <- data.frame(pine[[1]])
pine_tcontig <- data.frame(pine[[2]])

spruce <- species_data_import("all_spruce_data_pvalues.txt", TRUE)
overall_spruce <- data.frame(spruce[[1]])
spruce_tcontig <- data.frame(spruce[[2]])

cleaned_ortho <- ortholog_cleaning(pine_tcontig, spruce_tcontig, 3, oneToOne)


the_quantiles <- c (0.001,0.005,0.01,0.05)

binom_cuts <- c(0.99,0.999,0.9999,0.99999,0.999999,0.9999999,0.99999999,0.999999999)

pdf ("plot_top_candidate_tests.pdf")

for (i in 1:length (the_quantiles)){
  find_top_cand(pine, spruce, cleaned_ortho, the_quantiles, binom_cuts)
}


