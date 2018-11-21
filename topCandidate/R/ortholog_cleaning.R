#' ortholog_filtering
#'
#' Function to clean data for species entered
#' will get the unique names of value asked for, such as genes,in species then find which ones have above a certain number of occurances
#' Then returns a list of these values for the data provided 
#' @param species_data data to evaluate. For example, this could be the contigs present in a species
#' @param minimum the minimum number of occurences a value needs to have in this data to be included. For example, number of SNP's needed for a gene. Defaults to 3
#' @param columnName the column name which has the values to evaluate. For example, this could be the column of contigs. Defaults to NULL
#' @keywords 
#' @return returns the list of orthlogs that meet conditions and are in both species
#' @export
#' @examples
#' ortholog_cleaning(spruce_tcontig, pine_tcontig, 3)
#' 
#' pass name of column with the indicator variable as opposed to the tcontig, to make more general 
#' 

ortholog_filtering <- function(species_data, minimum = 3, columnName = NULL){
  columnLocation <- which( colnames(species_data)==columnName)
  #get the unique values of each of the column specified in input species
  values <- unique (as.character (species_one_tcontigs[[columnLocation]]))
  
  #with at least X occurences of test value:
  all <- table (as.character (species_one_tcontigs$[[columnLocation]]))
  numberFound <- all[all >= minimum]
  uniqueValues <- unique (names(numberFound))
  
  tabspr1 <- table (as.character (species_two_tcontigs$tcontig))
  tabspr2 <- tabspr1[tabspr1 >= minimum]
  spr_genes2 <- unique (names(tabspr2))
  
  
  #orthologs that are present in the data-tables of both species ##ADJUST as necessary depending on how many SNPs you want as a minimum
  int1 <- ortho[ortho$pine %in% pin_genes2,]
  ortho_in_both <- int1[int1$spruce %in% spr_genes2,]
  
  return(ortho_in_both)
}
