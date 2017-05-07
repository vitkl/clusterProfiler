##' GO term similarity - kappa score
##' Function will calculate kappa score between two categories describing a set of 
##' elements (GO terms, KEGG pathways, genes).
##' Function takes a data.table with 2 columns describing to two categories to 
##' which each element in rows can belong to. 
##' Os and 1s in the table signify membership.
##' defined as described in this article:
##' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2375021/figure/F2/
##' relies on data.table for subsetting and vcd package for calculating kappa 
##' score given 2x2 incidence matrix
##' the function is faster than irr::kappa2(ratings = value, weight = "unweighted")
##' the function is ~25 times slower than dist()
##'
##' @param value a data.table with 2 columns (two categories) and many rows (objects)
##' @importFrom vcd Kappa
##' @import data.table
##' @seealso \code{\link{categ_dist}}
##' @export
##' @author Vitalii Kleshchevnikov
kappa_score = function(value){
  colnames(value) = c("x","y")
  table_ = value[,.N, by = .(x,y)]
  table_2 = matrix(,2,2)
  table_2[1,1] = ifelse(length(table_[x == 1 & y == 1, N]), table_[x == 1 & y == 1, N],0)
  table_2[1,2] = ifelse(length(table_[x == 1 & y == 0, N]), table_[x == 1 & y == 0, N],0)
  table_2[2,1] = ifelse(length(table_[x == 0 & y == 1, N]), table_[x == 0 & y == 1, N],0)
  table_2[2,2] = ifelse(length(table_[x == 0 & y == 0, N]), table_[x == 0 & y == 0, N],0)
  kappa_score = vcd::Kappa(table_2,weights = "Fleiss-Cohen")$Unweighted
  return(kappa_score)
}

##' the categ_dist function to calculate categorical distance (Cohen's Kappa score) between multiple terms 
##' the function is intended to measure distances between GO terms based on proteins they annotate
##' more generally, the function can be used to measure categorical distances between any terms(categories) annotating objects
##' objects should be provided as a first column of a data.table, terms should be provided as a second column
##' 
##' Important: for correct evalutation mapping_table should inlude all 
##' term-to-gene associations, terms_to_compare should specify which you want to compare
##' @param mapping_table two-column data.table which provides mapping from object to term 
##' @param terms_to_compare a character vector specifying which terms to compare
##' @param ignore_limit logical, ignore the limit set on the number of terms
##' @import data.table
##' @importFrom caTools combs
##' @export
##' @author Vitalii Kleshchevnikov
categ_dist = function(mapping_table, terms_to_compare = unlist(unique(mapping_table[,2,with = F])), ignore_limit = F){
  if(ncol(mapping_table) > 2) stop("table has more than 2 columns, object id column and term column")
  if(ignore_limit == F) if(length(terms_to_compare) > 1000) stop("more than 1000 terms to compare, set ignore_limit = T if you are sure to proceed")
  if(!is.data.table(mapping_table)) stop("provided mapping / annotation table may not be in the right format (wrong class: not data.table)")
  #require(data.table)
  
  mapping_table = copy(unique(mapping_table))
  print(mapping_table)
  colnames(mapping_table) = c("UNIPROT", "GO")
  z1 = cbind(copy(mapping_table[,c("UNIPROT", "GO"), with = F]), value = 1)
  z2 = dcast(z1, UNIPROT ~ GO, fill = 0, drop = F)
  z2 = z2[,UNIPROT := NULL]
  z2 = z2[,terms_to_compare, with=F]
  combinations = t(combs(colnames(z2),2))
  dist = t(sapply(as.data.table(combinations), function(x) kappa_score(z2[,c(x[1],x[2]),with = F])))
  dist = cbind(as.data.table(dist), as.data.table(t(combinations)))
  colnames(dist) = c("kappa_score", "kappa_error", "GO1", "GO2")
  dist_temp = unique(rbind(dist,dist[,.(kappa_score,kappa_error, GO1 = GO2, GO2 = GO1)]))
  
  dist2 = as.matrix(dcast(dist_temp[,.(GO1,GO2, kappa_score)], GO1 ~ GO2))
  rownames_dist2 = dist2[,"GO1"]
  dist2 = as.matrix(dcast(dist_temp[,.(GO1,GO2, kappa_score)], GO1 ~ GO2)[,GO1 := NULL])
  rownames(dist2) = rownames_dist2
  dist2 = dist2[sort(rownames(dist2)), sort(colnames(dist2))]
  diag(dist2) = 1
  return(list(similarity_matrix = dist2, kappa_score_table = dist, kappa_score_table_redundant = dist_temp))
}