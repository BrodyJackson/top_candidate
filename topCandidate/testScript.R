

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

top_candidates <- find_top_candidate(pine, spruce, cleaned_ortho, the_quantiles, binom_cuts)

top_candidate_overlap <- top_candidates[[1]]
top_candidates_species_one <- top_candidates[[2]]
top_candidates_species_two <- top_candidates[[3]]


#calculate the hypergeometric statistics for one case (stephyper):

number_genes <- nrow (ortho_in_both) #total number of genes that can be considered as having been tested 
number_pine <- num_pin[3,3]
number_spruce <- num_spr[3,3]
number_overlap <- overlap[3,3]

#calculate the probability of getting number_overlap or more
sum(dhyper((number_overlap:number_pine),number_pine,(number_genes - number_pine),number_spruce))
