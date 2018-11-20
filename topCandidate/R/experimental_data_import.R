#' experimental_data_import
#'
#' Allows input of experimental info about a species
#' @param overall_info required table formatted data containing all the experimental info about a species. Defaults to NULL
#' @param get_only_annotated will get only the SNP's annotated to a gene if true, defaults to TRUE
#' @param annotatedColumn the column which you want to ensure that each value is annotated to. As an example this could be the contig column if you only wanted SNPS annotated to a gene. Defaults to NULL
#' @param commentChar the character which indicates comments in the input data frame, defaults to &
#' @keywords species
#' @export
#' @examples
#' experimental_data_import("all_pine_data_pvalues.txt", TRUE)
#' 
#' 

experimental_data_import <- function(overallInfo = NULL, get_only_annotated = TRUE, annotatedColumn = NULL, commentChar = "&"){
  #read in species data on p-values for each environmental test:
  if(is.null(overallInfo)){
    stop("Experimental info must be provided")
  }
  
  species_data <- read.table (overallInfo, T, comment.char = commentChar)
  
  for (i in 33:54){
    species_data[,i] <- as.numeric(as.character (species_data[,i]))
  }
  
  #get only the SNPs that are annotated to a gene (based on 3rd column)
  if(get_only_annotated == TRUE){
    if(is.null(annotatedColumn)){
      stop("Column to annotate values to must also be supplied if get_only_annotated == TRUE")
    }
    columnLocation <- which( colnames(species_data)==annotatedColumn)
    species_annotated <- species_data[grep("comp",species_data[[columnLocation]]),]
    returnList <- list(species_data, species_annotated)
  }
  else{
    returnList <- list(species_data)
  }
  ##assign the values read to the title entered as parameter
  return(returnList)

}

