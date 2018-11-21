#' find_top_candidates
#'
#' 
#' @param species_one_candidates list containing data frames holding top candidate values for various quantiles of species one. Defaults to NULL
#' @param species_two_candidates list containing data frames holding top candidate values for various quantiles of species two. Defaults to NULL
#' @param quantiles the quantiles being tested
#' @param binom the binomial cuts to test 
#' @param orthologs data frame containing the one to one orthologs of the value being tested for the two species 
#' @param ortho_in_both data frame containing the orthologs that are present in both the species being tests
#' @param species_one_name the name of species one (that is being used as column title)
#' @param species_two_name the name of species two (that is being used as column title)
#' @return 
#' @keywords compare
#' @export
#' @examples
#' find_top_candidates(pine, quantiles, tcontig, "MCMT")

compare_top_candidates <- function(species_one_candidates = NULL ,species_two_candidates = NULL, the_quantiles, binom, orthologs, ortho_in_both, species_one_name = "pine", species_two_name = "spruce"){
  
  if(is.null(species_one_candidates) || is.null(species_two_candidates)){
    stop("need top candidate data of two species")
  }
 
  

  overlap <- array (NA, c (length (the_quantiles),length (binom)))
  num_species_one <- array (NA, c (length (the_quantiles),length (binom)))
  num_species_two <- array (NA, c (length (the_quantiles),length (binom)))
  num_ortho <- array (NA, c (length (the_quantiles),length (binom)))
  
  
  pdf ("top_candidates_compared.pdf")
  
  
  
  for (i in 1:length (the_quantiles)){
    
    
    sub_good_species_one <- species_one_candidates[[i]]
    sub_good_species_one <- merge (sub_good_species_one,orthologs,by.x = "names.snps_count2.",by.y = species_one_name, all.x = T) #many of the rows do not have an ortholog identified in the other species, so there is an "NA" in the 4th column
    
    totsnp1 <- sum (sub_good_species_one$snps_count2)
    totout1 <- sum (sub_good_species_one$outliers_count2)
    expect1_species_one <- totout1 / totsnp1 #expected number of outliers per SNP, based on the total number of outliers and the total number of SNPs in the genes that have at least one outlier. Other choices are possible here.
    
    
    sub_good_species_two <- species_two_candidates[[i]]
    sub_good_species_two <- merge (sub_good_species_two,orthologs,by.x = "names.snps_count2.",by.y = species_two_name, all.x = T) #many of the rows do not have an ortholog identified in the other species, so there is an "NA" in the 4th column
    
    totsnp1 <- sum (sub_good_species_two$snps_count2)
    totout1 <- sum (sub_good_species_two$outliers_count2)
    expect1_species_two <- totout1 / totsnp1 #expected number of outliers per SNP, based on the total number of outliers and the total number of SNPs in the genes that have at least one outlier. Other choices are possible here.
    
    
    
    for (j in 1:length (binom_cuts)){
      
      
      #calculate all of the top candidate test cutoffs
      binom_species_one <- qbinom (binom_cuts[j],sub_good_species_one$snps_count2, expect1_species_one)	
      binom_species_two <- qbinom (binom_cuts[j],sub_good_species_two$snps_count2, expect1_species_two)	
      
      #identify which cases have more outliers than the cutoff
      sub_species_one <- sub_good_species_one[which (sub_good_species_one$outliers_count2 > binom_species_one),]
      sub_species_two <- sub_good_species_two[which (sub_good_species_two$outliers_count2 > binom_species_two),]
      
      #Could be error here, need to change the spruce and pine selectors to more general
      sub_species_one <- sub_species_one[which (is.na(sub_species_one$spruce) == F),]
      sub_species_two <- sub_species_two[which (is.na(sub_species_two$pine) == F),]
      
      column_location1 <- which( colnames(ortho_in_both)==species_one_name)
      column_location2 <- which( colnames(ortho_in_both)==species_two_name)
      column_location_sub <- which( colnames(sub_species_one)==species_two_name)
      
      overlap[i,j] <- sum (sub_species_one[[column_location_sub]] %in% sub_species_two$names.snps_count2.) #how many genes are top candidates in both species?
      num_species_one[i,j] <- sum (sub_species_one$names.snps_count2. %in% ortho_in_both[[column_location1]]) #how many genes are top candidates in pine?
      num_species_two[i,j] <- sum (sub_species_two$names.snps_count2. %in% ortho_in_both[[column_location2]]) #how many genes are top candidates in spruce?
      
      par (mfcol = c (2,1))
      par (mar = c (5,5,5,5))
      plot (sub_good_pine$snps_count2,sub_good_pine$outliers_count2, xlab = "number of SNPs", ylab = "number of outliers", main = paste ("quantile = ", the_quantiles[i],"  binom_cut = ",binom[j]),  cex = 0.5)
      points (1:(max(sub_good_species_one$snps_count2)+10),qbinom(binom[j],1:(max(sub_good_species_one$snps_count2)+10), expect1_species_one), type = "l", col = "red")
      
      #merge together the pine and spruce tables:
      sub_good_species_one_2 <- cbind (sub_good_species_one,binom_species_one)
      sub_good_species_two_2 <- cbind (sub_good_species_two,binom_species_two)
      
      ##values not changed to species1/2 properly yet
      # merg1 <- merge (sub_good_species_one_2,sub_good_species_two_2,by.x = "names.snps_count2.",by.y = "pine")
      # par (mar = c (5,5,5,5))
      # plot ((merg1$outliers_count2.x / merg1$binom_pi),(merg1$outliers_count2.y / merg1$binom_sp), xlab = "Pine",ylab = "Spruce", cex = 0.5,xlim = c (0,5),ylim = c (0,5), main = "Number of outliers / cutoff")
      # arrows (-1000,1,1000,1)
      # arrows (1,-1000,1,1000)
    }	
    
  }
  dev.off()
  
  number_genes <- nrow (ortho_in_both) #total number of genes that can be considered as having been tested 
  number_species_one <- num_species_one[3,3]
  number_species_two <- num_species_two[3,3]
  number_overlap <- overlap[3,3]
  
  #calculate the probability of getting number_overlap or more
  sum(dhyper((number_overlap:number_pine),number_pine,(number_genes - number_pine),number_spruce))
}