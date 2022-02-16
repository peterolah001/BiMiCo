
#' Removal of non-bacterial sequences from a phyloseq object
#'
#' Given an input phyloseq object with valid tax_table, returns taxonomy-filtered phyloseq object with mitochondrial, chloroplast and non-bacterial taxa removed.
#' Non-bacterial sequences (e.g. mitochondria, chloroplast) need to be removed from OTU/ASV tables prior to statistical analyses. rm_nonbac removes all rows of otu_table and tax_table in a phyloseq object based on taxonomic labelling 'Mitochondria', 'Chloroplast' or where Kingdom is other than 'Bacteria'. NAs in taxonomy table are converted to 'unidentified'.  
#' @import phyloseq
#' @param x Phyloseq object with non-empty otu_table and tax_table slots; otu_table may contain either OTUs or ASVs. 
#' @keywords data processing phyloseq
#' @export
#' @examples
#' rm_nonbac()

rm_nonbac = function(x){
  if (any(is.na(tax_table(x)))){
    tax_table(x)[is.na(tax_table(x))]<-"unidentified"
  } else {
    tax_table(x) <- x
  }
  if (any(tax_table(x)[,1]!="Bacteria" )) {
    nonbac = c(taxa_names(subset_taxa(x, Kingdom!="Bacteria")))
  } else {
    nonbac = as.character()
  }
  if (any(tax_table(x)[,3]=="Chloroplast", na.rm=TRUE)) {
    chl = c(taxa_names(subset_taxa(x, Class=="Chloroplast")))
  }else {
    chl = as.character()
  }
  if (any(tax_table(x)[,6]=="Mitochondria", na.rm=TRUE)) {
    mito = c(taxa_names(subset_taxa(x, Family=="Mitochondria")))
  }else {
    mito = as.character()
  }
  removetaxa = c(nonbac, chl, mito)
  allTaxa = taxa_names(x)
  myTaxa <- allTaxa[!(allTaxa %in% removetaxa)]
  new_phy = prune_taxa(myTaxa, x)
  return(new_phy)
}
