#' ortholog_cleaning
#'
#' Function to clean SNP data for two species entered
#' will get the unique names of each gene in species then find which ones have above a certain number of SNP's
#' Then returns the list representing the ortologs meeting this criteria that are in both species
#' @param species_one_tcontigs contig data from species one
#' @param species_two_tcontigs contig data from species two
#' @param minimum the minimum number of SNP's needed for a gene
#' @param one_to_one_orthologs Data frame representing the one to one orthologs of the two species
#' @keywords clean
#' @return returns the list of orthlogs that meet conditions and are in both species
#' @export
#' @examples
#' ortholog_cleaning(spruce_tcontig, pine_tcontig, 3)

ortholog_cleaning <- function(species_one_tcontigs, species_two_tcontigs, minimum, ortho){
  #get the unique names of each of the genes in each species
  pin_genes1 <- unique (as.character (species_one_tcontigs$tcontig))
  spr_genes1 <- unique (as.character (species_two_tcontigs$tcontig))
  
  #with at least X SNPs:
  tabpin1 <- table (as.character (species_one_tcontigs$tcontig))
  tabpin2 <- tabpin1[tabpin1 >= minimum]
  pin_genes2 <- unique (names(tabpin2))
  
  tabspr1 <- table (as.character (species_two_tcontigs$tcontig))
  tabspr2 <- tabspr1[tabspr1 >= minimum]
  spr_genes2 <- unique (names(tabspr2))
  
  
  #orthologs that are present in the data-tables of both species ##ADJUST as necessary depending on how many SNPs you want as a minimum
  int1 <- ortho[ortho$pine %in% pin_genes2,]
  ortho_in_both <- int1[int1$spruce %in% spr_genes2,]
  
  return(ortho_in_both)
}
