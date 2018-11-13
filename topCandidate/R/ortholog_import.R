#' ortholog_import
#'
#' This function allows you to import all the necessary ortholog information between two species
#' Takes two table of ortholog sequences, one of which is all orthologs, and the other just the one to one
#' @param all_orthologs table formatted data containing all the orthologs between two species
#' @param one_to_one_orthologs table formatted data containine the one to one orthologs of these species
#' @keywords orthologs
#' @return returns a list of data Frames of the read in ortholog data
#' @export
#' @examples
#' ortholog_import("all_orthologs_pine_spruce.txt","one_to_one_orthologs.txt")

ortholog_import <- function(all_orthologs, one_to_one_orthologs){
  all_ortho <- read.table (all_orthologs,T)
  ortho <- read.table (one_to_one_orthologs, T)
  returnList <- list(all_ortho, ortho)
  return(returnList)
}

