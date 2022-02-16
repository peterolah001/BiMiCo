#' Convert Phyloseq object with ASVs renamed to more convenient shorter form, and actual ASV sequences stored as a DNAStringSet (refseq slot)
#'
#' Takes as input a phyloseq object and converts otu_table and tax_table rownames to shorter form, stores ASV sequences in refseq slot as 'DNAStringSet'
#' @import phyloseq
#' @param phylo (Required) input object of class 'phyloseq'
#' @keywords read processing dada2
#' @export
#' @examples
#' rename_asvs()

rename_asvs <- function(phylo){
  
  seqs <- Biostrings::DNAStringSet(taxa_names(phylo))
  names(seqs) <- phyloseq::taxa_names(phylo)
  phs <- phyloseq::merge_phyloseq(phylo, seqs)
  phyloseq::taxa_names(phylo) <- paste0("ASV", seq(phyloseq::ntaxa(phylo)))
  
  return(phylo)
  
}

