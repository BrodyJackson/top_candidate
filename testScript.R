
##Script showiing example package usage
##seperate top candidate into each species, then allow a compare function 
##specify outliers as upper or lower
## fst vs p value 

install("topCandidate")
library(topCandidate)


##read in correspondence between pine and spruce orthologs, including both one-to-one and multi-to one correspondences
allOrthologs <- read.table ("all_orthologs_pine_spruce.txt",T)
## read in one-to-one orthologs only:
oneToOne  <- read.table ("one_to_one_orthologs.txt", T)

pine <- experimental_data_import("all_pine_data_pvalues.txt", TRUE, "tcontig", "&")
overall_pine <- data.frame(pine[[1]])
pine_tcontig <- data.frame(pine[[2]])

spruce <- experimental_data_import("all_spruce_data_pvalues.txt", TRUE, "tcontig", "&")
overall_spruce <- data.frame(spruce[[1]])
spruce_tcontig <- data.frame(spruce[[2]])

filtered_orthologs_pine <- ortholog_filtering(pine_tcontig, 3, "tcontig")
filtered_orthologs_spruce <- ortholog_filtering(spruce_tcontig, 3, "tcontig")

#orthologs that are present in the data-tables of both species ##ADJUST as necessary depending on how many SNPs you want as a minimum
#This part can be added when the compare function is created
##int1 <- ortho[ortho$pine %in% pin_genes2,]
##ortho_in_both <- int1[int1$spruce %in% spr_genes2,]

##allows changing of quantiles and binomial cuts to be used when finding top candidates for a species
the_quantiles <- c (0.001,0.005,0.01,0.05)

binom_cuts <- c(0.99,0.999,0.9999,0.99999,0.999999,0.9999999,0.99999999,0.999999999)

#finds the outliers for each species and returns a list, list[1] contains number of SNPs with at least one outlier
#list[2] contains number of outliers
#results can be passed into top candidate test to compare
pine_outliers <- find_top_candidates(pine, the_quantiles, binom_cuts, "tcontig", "MCMT")

spruce_outliers <- find_top_candidates(spruce, the_quantiles, binom_cuts, "tcontig", "MCMT")



# find_top_candidate(pine_outliers, spruce_outliers, binom_cuts, the_quantiles)
# # 
# # top_candidate_overlap <- top_candidates[[1]]
# # top_candidates_species_one <- top_candidates[[2]]
# # top_candidates_species_two <- top_candidates[[3]]
# 
# 
# #calculate the hypergeometric statistics for one case (stephyper):
# 
# number_genes <- nrow (ortho_in_both) #total number of genes that can be considered as having been tested 
# number_pine <- num_pin[3,3]
# number_spruce <- num_spr[3,3]
# number_overlap <- overlap[3,3]
# 
# #calculate the probability of getting number_overlap or more
# sum(dhyper((number_overlap:number_pine),number_pine,(number_genes - number_pine),number_spruce))
