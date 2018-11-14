#' find_top_candidate
#'
#' 
#' @param species_one_data list containing overall data in first position and tcontigs in second, for first species
#' @param species_two_data list containing overall data in first position and tcontigs in second, for second species
#' @param orthologs orthologs present in both, which have been cleaned
#' @param quantiles the quantiles being tested
#' @param binomial_cuts the binomial cuts to test
#' @return returns a list containing the top candidates. First position is candidate overlap, second is candidates in species one, third is candidates in species 2
#' @keywords topcandidate
#' @export
#' @examples
#' find_top_cand(pine, spruce, cleaned_ortho, binomial_cuts)

find_top_candidate <- function(species_one_data, species_two_data, orthologs, the_quantiles, binomial_cuts){
  
  pdf ("plot_top_candidate_tests.pdf")
  
  overlap <- array (NA, c (length (the_quantiles),length (binomial_cuts)))
  num_species_one <- array (NA, c (length (the_quantiles),length (binomial_cuts)))
  num_species_two <- array (NA, c (length (the_quantiles),length (binomial_cuts)))
  
  for (i in 1:length (the_quantiles)){
    
    the_quantile <- the_quantiles[i]
    
    #PINE: calculate the quantile to call outliers over all environmental variables. This somewhat arbitrary choice allows some environmental variables to have more outliers than others, based on the relative strength of the p-values
    quantile1 <- quantile (species_one_data[[1]][,33:54], the_quantile,na.rm = T)	#STEP2
    
    
    the_i <- which (colnames (species_one_data[[2]]) == "MCMT") 
    
    outliers <- species_one_data[[2]][,the_i] < quantile1
    snps <- species_one_data[[2]][,the_i] < 1000000  ## arbitrary, just allows using the same code as outliers
    
    
    outliers_count <- tapply (outliers, list(as.character (species_one_data[[2]]$tcontig)),sum, na.rm = T)  #count the number of outliers per gene
    outliers_count2 <- outliers_count[outliers_count >=1] ## only those that have at least one outlier
    
    snps_count <- tapply (snps, list(as.character (species_one_data[[2]]$tcontig)),sum, na.rm = T) #count the number of SNPs per gene
    snps_count2 <- snps_count[outliers_count >=1] ## only those that have at least one outlier
    
    sub_good_pine <- data.frame (names(snps_count2),snps_count2,outliers_count2)
    sub_good_pine <- merge (sub_good_pine,orthologs,by.x = "names.snps_count2.",by.y = "pine", all.x = T) #many of the rows do not have an ortholog identified in the other species, so there is an "NA" in the 4th column
    
    
    totsnp1 <- sum (sub_good_pine$snps_count2)
    totout1 <- sum (sub_good_pine$outliers_count2)
    expect1_pi <- totout1 / totsnp1 #expected number of outliers per SNP, based on the total number of outliers and the total number of SNPs in the genes that have at least one outlier. Other choices are possible here.
    
    #NOTE: the outlier threshold expect1_pi is much greater than 1%. This is because the quantile is calculated over all environmental variables. 
    
    #SPRUCE:
    quantile1 <- quantile (species_two_data[[1]][,33:54], the_quantile,na.rm = T)		#STEP2
    
    the_i <- which (colnames (species_two_data[[2]]) == "MCMT") 
    
    outliers <- species_two_data[[2]][,the_i] < quantile1
    snps <- species_two_data[[2]][,the_i] < 1000000
    
    outliers_count <- tapply (outliers, list(as.character (species_two_data[[2]]$tcontig)),sum, na.rm = T) #count the number of outliers per gene
    outliers_count2 <- outliers_count[outliers_count >=1]
    
    snps_count <- tapply (snps, list(as.character (species_two_data[[2]]$tcontig)),sum, na.rm = T) #count the number of SNPs per gene
    snps_count2 <- snps_count[outliers_count >=1]
    
    sub_good_spruce <- data.frame (names(snps_count2),snps_count2,outliers_count2)
    sub_good_spruce <- merge (sub_good_spruce,orthologs,by.x = "names.snps_count2.",by.y = "spruce", all.x = T) #many of the rows do not have an ortholog identified in the other species, so there is an "NA" in the 4th column
    
    
    totsnp1 <- sum (sub_good_spruce$snps_count2)  
    totout1 <- sum (sub_good_spruce$outliers_count2)
    expect1_sp <- totout1 / totsnp1
    
    for (j in 1:length (binomial_cuts)){
      
      #calculate all of the top candidate test cutoffs
      binom_pi <- qbinom (binomial_cuts[j],sub_good_pine$snps_count2, expect1_pi)	
      binom_sp <- qbinom (binomial_cuts[j],sub_good_spruce$snps_count2, expect1_sp)	
      
      #identify which cases have more outliers than the cutoff
      sub_pin <- sub_good_pine[which (sub_good_pine$outliers_count2 > binom_pi),]
      sub_spr <- sub_good_spruce[which (sub_good_spruce$outliers_count2 > binom_sp),]
      
      sub_pin <- sub_pin[which (is.na(sub_pin$spruce) == F),]
      sub_spr <- sub_spr[which (is.na(sub_spr$pine) == F),]
      
      overlap[i,j] <- sum (sub_pin$spruce %in% sub_spr$names.snps_count2.) #how many genes are top candidates in both species?
      num_species_one[i,j] <- sum (sub_pin$names.snps_count2. %in% orthologs$pine) #how many genes are top candidates in pine?
      num_species_two[i,j] <- sum (sub_spr$names.snps_count2. %in% orthologs$spruce) #how many genes are top candidates in spruce?
      
      par (mfcol = c (2,1))
      par (mar = c (5,5,5,5))
      plot (sub_good_pine$snps_count2,sub_good_pine$outliers_count2, xlab = "number of SNPs", ylab = "number of outliers", main = paste ("quantile = ", the_quantiles[i],"  binom_cut = ",binomial_cuts[j]),  cex = 0.5)
      points (1:(max(sub_good_pine$snps_count2)+10),qbinom(binomial_cuts[j],1:(max(sub_good_pine$snps_count2)+10), expect1_pi), type = "l", col = "red")
      
      #merge together the pine and spruce tables:
      sub_good_pine2 <- cbind (sub_good_pine,binom_pi)
      sub_good_spruce2 <- cbind (sub_good_spruce,binom_sp)
      
      merg1 <- merge (sub_good_pine2,sub_good_spruce2,by.x = "names.snps_count2.",by.y = "pine")
      par (mar = c (5,5,5,5))
      plot ((merg1$outliers_count2.x / merg1$binom_pi),(merg1$outliers_count2.y / merg1$binom_sp), xlab = "Pine",ylab = "Spruce", cex = 0.5,xlim = c (0,5),ylim = c (0,5), main = "Number of outliers / cutoff")
      arrows (-1000,1,1000,1)
      arrows (1,-1000,1,1000)
    }	
  }
  dev.off()
  
  returnList <- list(overlap, num_species_one, num_species_two)
  return(returnList)
}