install("topCandidateTest")
library(topCandidateTest)


##read in correspondence between pine and spruce orthologs, including both one-to-one and multi-to one correspondences
allOrthologs <- read.table ("all_orthologs_pine_spruce.txt",T)
## read in one-to-one orthologs only:
oneToOne  <- read.table ("one_to_one_orthologs.txt", T)

pine <- experimental_data_import("all_pine_data_pvalues.txt", TRUE, "tcontig")
overall_pine <- data.frame(pine[[1]])
pine_tcontig <- data.frame(pine[[2]])

spruce <- experimental_data_import("all_spruce_data_pvalues.txt", TRUE, "tcontig")
overall_spruce <- data.frame(spruce[[1]])
spruce_tcontig <- data.frame(spruce[[2]])

filtered_orthologs_pine <- ortholog_filtering(pine_tcontig, 3, "tcontig")
filtered_orthologs_spruce <- ortholog_filtering(spruce_tcontig, 3, "tcontig")

##Could try and implement this part as a seperate function? 
int1 <- oneToOne[oneToOne$pine %in% filtered_orthologs_pine,]
ortho_in_both <- int1[int1$spruce %in% filtered_orthologs_spruce,]

##allows changing of quantiles and binomial cuts to be used when finding top candidates for a species
the_quantiles <- c (0.001,0.005,0.01,0.05)

binom_cuts <- c(0.99,0.999,0.9999,0.99999,0.999999,0.9999999,0.99999999,0.999999999)

pine_outliers <- find_top_candidates(pine, the_quantiles, binom_cuts, "tcontig", "MCMT", "low", "pine_top_candidates.pdf")

spruce_outliers <- find_top_candidates(spruce, the_quantiles, binom_cuts, "tcontig", "MCMT", "low", "spruce_top_candidates.pdf")

# compare top candidates will return a matrix holding the probability of getting top candidates overlapping in both species for each quantile and binomial cut comparison
# the quantiles are rows while the binomial cuts are columns 
probability <- compare_top_candidates(pine_outliers, spruce_outliers, the_quantiles, binom_cuts, oneToOne, ortho_in_both, "pine", "spruce", "top_candidates_compared.pdf")
#the probability of getting a certain number of top candidates overlapping in both species
print(probability)

