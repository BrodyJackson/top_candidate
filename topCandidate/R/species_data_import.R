#' species_data_import
#'
#' Allows input of experimental info about a species
#' @param overall_info table formatted data containing all the experimental info about a species
#' @param get_only_annotated will get only the SNP's annotated to a gene if true, defaults to TRUE
#' @keywords species
#' @export
#' @examples
#' species_data_import("all_pine_data_pvalues.txt", TRUE)

species_data_import <- function(overallInfo, get_only_annotated = TRUE){
  #read in species data on p-values for each environmental test:
  species_data <- read.table (overallInfo, T, comment.char = "&")
  
  for (i in 33:54){
    species_data[,i] <- as.numeric(as.character (species_data[,i]))
  }
  
  #get only the SNPs that are annotated to a gene (based on 3rd column)
  if(get_only_annotated == TRUE){
    species_tcontig <- species_data[grep("comp",species_data$tcontig),]
    returnList <- list(species_data, species_tcontig)
  }
  else{
    returnList <- list(species_data)
  }
  ##assign the values read to the title entered as parameter
  return(returnList)

}

