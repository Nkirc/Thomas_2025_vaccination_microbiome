#############################
# Taxonomy string splitter #
############################

#come back to string splitter and make it shorter (function in a function?)
#function for each individual unsplit group - and can then have a vector of QIIME notation or something
string_splitter  <- function(tax_col){
  #takes a one column matrix of qiime taxonomy output and splits individual taxonomic levels into their own columns
  strain <- tax_col
  species <- tax_col
  genus <- tax_col
  family <- tax_col
  order <- tax_col
  class <- tax_col
  phylum <- tax_col
  kingdom <- tax_col
  for(i in 1:length(tax_col)){
    #Split before domain name
    kingdom[i]     <- unlist(strsplit(kingdom[i], split='k__',   fixed=T))[2]   #Split after 'k__' in taxa string
    #Split after domain name
    kingdom[i]     <- unlist(strsplit(kingdom[i], split='|p__', fixed=T))[1]   #Split before '|p__' in taxa string
    #Split before phylum name
    phylum[i]     <- unlist(strsplit(phylum[i], split='p__',   fixed=T))[2]   #Split after 'p__' in taxa string
    #Split after phylum name
    phylum[i]     <- unlist(strsplit(phylum[i], split='|c__', fixed=T))[1]   #Split before '|c__' in taxa string
    #Split before class name
    class[i]     <- unlist(strsplit(class[i], split='c__',   fixed=T))[2]   #Split after 'c__' in taxa string
    #Split after class name
    class[i]     <- unlist(strsplit(class[i], split='|o__', fixed=T))[1]   #Split before '|f__' in taxa string
    #Split before class name
    order[i]     <- unlist(strsplit(order[i], split='o__',   fixed=T))[2]   #Split after 'o__' in taxa string
    #Split after class name
    order[i]     <- unlist(strsplit(order[i], split='|f__', fixed=T))[1]   #Split before '|f__' in taxa string
    #Split before family name
    family[i]     <- unlist(strsplit(family[i], split='f__',   fixed=T))[2]   #Split after 'f__' in taxa string
    #Split after family name
    family[i]     <- unlist(strsplit(family[i], split='|g__', fixed=T))[1]   #Split before '|g__' in taxa string
    #Split before genus name
    genus[i]     <- unlist(strsplit(genus[i], split='g__',   fixed=T))[2]   #Split after 'g__' in taxa string
    #Split after genus nam|
    genus[i]     <- unlist(strsplit(genus[i], split='|s__', fixed=T))[1]   #Split before '|s__' in taxa string
    #Split before species name
    species[i]     <- unlist(strsplit(species[i], split='s__',   fixed=T))[2]   #Split after 's__' in taxa string
    #Split after species name
    species[i]     <- unlist(strsplit(species[i], split='|t__', fixed=T))[1]
    #Split before strain name
    strain[i]     <- unlist(strsplit(strain[i], split='t__',   fixed=T))[2]
  }
  all_tax <- cbind(kingdom, phylum, class, order, family, genus, species, strain)
  colnames(all_tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  return(all_tax)
}