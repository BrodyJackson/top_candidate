#' find_top_candidates
#'
#' 
#' @param species_data list containing overall data in first position and orthologs data to evaluate in second. These values are previously found using experimental_data_import method. Defaults to NULL
#' @param quantiles the quantiles being tested
#' @param binom the binomial cuts to test 
#' @param candidateColumnName the name of the column containing all the candidates that will be tested. Defaults to NULL
#' @param testVariable the column name holding the results for a certain test. For example, this could be the results for a certain environmental location
#' @param highOrLow value to indicate if experimental results should be higher or lower than the quantile value being tested. For example, would be changed depending on if results are p-value or FST. Defaults to low
#' @param plotFileName the name of the pdf file you want to plot in, if NULL then top candidates will not be plotted, defaults to null
#' @return returns a list containing data frames holding the top candidate results for each quantile 
#' @keywords topcandidate
#' @export
#' @examples
#' find_top_candidates(pine, the_quantiles, binom_cuts, "tcontig", "MCMT", "pine_top_candidates.pdf")

find_top_candidates <- function(species_data = NULL ,the_quantiles, binom, candidateColumnName = NULL, testVariable = NULL, highOrLow = "low", plotFileName = NULL){
  
  if(is.null(candidateColumnName)){
    stop("Column name of candidates must be supplied")
  }
  if(is.null(species_data)){
    stop("species data must be provided")
  }
  if(is.null(testVariable)){
    stop("test variable must be provided")
  }
  if((highOrLow != "low") && (highOrLow != "high")){
    stop("Value of highOrLow must be supplied as either 'high' or 'low' for testing experimental condition")
  }
  
  returnList <- list()
  
  columnLocation <- which( colnames(species_data)==candidateColumnName)
  
  num_found <- array (NA, c (length (the_quantiles),length (binom)))
  
  if(!is.null(plotFileName)){
    pdf (plotFileName)
  }
  
  for (i in 1:length (the_quantiles)){
    
    the_quantile <- the_quantiles[i]
    
    the_i <- which (colnames (species_data[[2]]) == testVariable)
    
    #calculate the quantile to call outliers over all environmental variables. Used only on the testColumn name supplied
    quantile1 <- quantile (species_data[[1]][[the_i]], the_quantile,na.rm = T)	#STEP2
    
    # specify weather we are looking for values above or beloe the quantile 
    if(highOrLow == "low"){
      outliers <- species_data[[2]][,the_i] < quantile1
      snps <- species_data[[2]][,the_i] < 1000000  ## arbitrary, just allows using the same code as outliers  
    }
    else if(highOrLow == "high"){
      outliers <- species_data[[2]][,the_i] > quantile1
      snps <- species_data[[2]][,the_i] < 1000000  
    }
    outliers_count <- tapply (outliers, list(as.character (species_data[[2]]$tcontig)),sum, na.rm = T)  #count the number of outliers per gene
    outliers_count2 <- outliers_count[outliers_count >=1] ## only those that have at least one outlier
    
    snps_count <- tapply (snps, list(as.character (species_data[[2]]$tcontig)),sum, na.rm = T) #count the number of SNPs per gene
    snps_count2 <- snps_count[outliers_count >=1] ## only those that have at least one outlier
    
    sub_good <- data.frame (names(snps_count2),snps_count2,outliers_count2)
    # sub_good_spruce <- merge (sub_good_spruce,ortho,by.x = "names.snps_count2.",by.y = "spruce", all.x = T) #many of the rows do not have an ortholog identified in the other species, so there is an "NA" in the 4th column
    returnList[[i]] <- sub_good
    
    totsnp1 <- sum (sub_good$snps_count2)
    totout1 <- sum (sub_good$outliers_count2)
    expected <- totout1 / totsnp1 #expected number of outliers per SNP, based on the total number of outliers and the total number of SNPs in the genes that have at least one outlier. Other choices are possible here.
    
    if(!is.null(plotFileName)){
      for (j in 1:length (binom)){
        par (mfcol = c (2,1))
        par (mar = c (5,5,5,5))
        plot (sub_good$snps_count2,sub_good$outliers_count2, xlab = "number of SNPs", ylab = "number of outliers", main = paste ("quantile = ", the_quantiles[i],"  binom_cut = ",binom[j]),  cex = 0.5)
        points (1:(max(sub_good$snps_count2)+10),qbinom(binom[j],1:(max(sub_good$snps_count2)+10), expected), type = "l", col = "red")
      }	
    }
  }
  dev.off()
  return(returnList)
}