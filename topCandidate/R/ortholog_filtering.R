#' ortholog_filtering
#'
#' Function to clean data for species entered
#' will get the unique names of value asked for, such as genes,in species then find which ones have above a certain number of occurances
#' Then returns a list of these values for the data provided 
#' @param species_data data to evaluate. For example, this could be the contigs present in a species. Defaults to NULL
#' @param minimum the minimum number of occurences a value needs to have in this data to be included. For example, number of SNP's needed for a gene. Defaults to 3
#' @param columnName the column name which has the values to evaluate. For example, this could be the column of contigs. Defaults to NULL
#' @keywords 
#' @return returns the list of orthlogs that meet conditions 
#' @export
#' @examples
#' ortholog_filtering(pine_tcontig, 3, "tcontig")
#' 
#' pass name of column with the indicator variable as opposed to the tcontig, to make more general 
#' 

ortholog_filtering <- function(species_data = NULL, minimum = 3, columnName = NULL){
  
  if(is.null(columnName)){
    stop("Column of values to examine must be supplied")
  }
  if(is.null(speciesData)){
    stop("Data to examine must be supplied")
  }
  
  columnLocation <- which( colnames(species_data)==columnName)
  #get the unique values of each of the column specified in input species
  values <- unique (as.character (species_one_tcontigs[[columnLocation]]))
  
  #with at least X occurences of test value:
  all <- table (as.character (species_one_tcontigs[[columnLocation]]))
  numberFound <- all[all >= minimum]
  filteredValues <- unique (names(numberFound))
  
  return(filteredValues)
}
