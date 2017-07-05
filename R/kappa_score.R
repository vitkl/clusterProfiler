##' similarity between categories - kappa score
##' 
##' @description \code{kappa_score} calculates kappa similarity score between two categories based on the overlap of the set of elements they describe adjusted for that overlap being produced by chance. 
##' @details Categories can be anything, common example is GO terms (describe gene products, \url{http://geneontology.org/}), KEGG pathways (describe genes, \url{http://www.genome.jp/kegg/pathway.html}), genes (similarity based on annotation with GO terms).
##' @details \code{kappa_score} takes a data.table with 2 columns describing to two categories to which each element in rows can belong to. Os and 1s in the table signify membership.
##' @details \code{kappa_score} is defined as described in this article: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2375021/figure/F2/}
##' @details relies on data.table for subsetting and vcd package for calculating kappa score given 2x2 incidence matrix
##' 
##' @details \code{kappa_score} is faster than irr::kappa2(ratings = value, weight = "unweighted")
##' @details \code{kappa_score} is ~25 times slower than dist()
##'
##' @param value a data.table with 2 columns (two categories) and many rows (objects)
##' @importFrom vcd Kappa
##' @import data.table
##' @seealso \code{\link{categ_dist}}
##' @export
##' @author Vitalii Kleshchevnikov
kappa_score = function(value){
    if(!is.data.table(value)) stop("kappa_score: provided table may not be in the right format (wrong class: not data.table)")
    if(ncol(value) > 2) stop("kappa_score: table has more than 2 category columns")
    
    colnames(value) = c("x","y")
    table_ = value[,.N, by = .(x,y)]
    table_2 = matrix(,2,2)
    # calculate the number of cases when the item belongs to both categories
    table_2[1,1] = ifelse(length(table_[x == 1 & y == 1, N]), table_[x == 1 & y == 1, N],0)
    # calculate the number of cases when the item belongs to x but not y
    table_2[1,2] = ifelse(length(table_[x == 1 & y == 0, N]), table_[x == 1 & y == 0, N],0)
    # calculate the number of cases when the item belongs to y but not x
    table_2[2,1] = ifelse(length(table_[x == 0 & y == 1, N]), table_[x == 0 & y == 1, N],0)
    # calculate the number of cases when the item belongs to neither x nor y
    table_2[2,2] = ifelse(length(table_[x == 0 & y == 0, N]), table_[x == 0 & y == 0, N],0)
    # calculate kappa score using Kappa function from vcd package
    kappa_score = Kappa(table_2,weights = "Fleiss-Cohen")$Unweighted
    return(kappa_score)
}

##' categorical distance (Cohen's Kappa score) between multiple categories
##' 
##' @description The \code{categ_dist} function calculates categorical distance (Cohen's Kappa score) between multiple categories (terms).
##' @details \code{categ_dist} is intended to measure distances between GO terms based on proteins they annotate.
##' @details More generally, \code{categ_dist} can be used to measure categorical distances between any terms(categories) annotating objects.
##' @details Objects should be provided as a first column of a data.table, terms (categories) should be provided as a second column.
##' 
##' @details Important: for correct evalutation of GO term similarity mapping_table should inlude all term-to-gene associations (expanded so that parent term inherits all annotations of children terms), terms_to_compare should specify which you want to compare
##' @details For full gene to GO mapping table (expanded so that parent term inherits all annotations of children terms) blows up the memory (over 50GB)
##' @param mapping_table two-column data.table which provides mapping from object to term 
##' @param terms_to_compare a character vector specifying which terms to compare
##' @param ignore_limit logical, ignore the limit set on the number of terms to compare (1000)
##' @param parallel logical, if TRUE use parSapply (from the parallel package) to speed up calculations on multiple cores
##' @import data.table
##' @importFrom caTools combs
##' @importFrom parallel detectCores
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterEvalQ
##' @importFrom parallel parSapply
##' @importFrom parallel stopCluster
##' @seealso \code{\link{kappa_score}}
##' @seealso \code{\link{parSapply}}
##' @export
##' @author Vitalii Kleshchevnikov
categ_dist = function(mapping_table, terms_to_compare = unlist(unique(mapping_table[,2,with = F])), ignore_limit = F, parallel = T){
    if(ncol(mapping_table) > 2) stop("categ_dist: table has more than 2 columns, object id column and term column")
    if(ignore_limit == F) if(length(terms_to_compare) > 1000) stop("categ_dist: more than 1000 terms to compare, set ignore_limit = T if you are sure to proceed")
    if(!is.data.table(mapping_table)) stop("categ_dist: provided mapping / annotation table may not be in the right format (wrong class: not data.table)")
    
    # copy original data to avoid allowing data.table to modify original data
    mapping_table = copy(unique(mapping_table))
    # sanity check: print(mapping_table)
    # assign standard column names
    colnames(mapping_table) = c("UNIPROT", "GO")
    # covert two column table into N UNIPROT by N GO table containing 1 when UNIPROT-1 is in group GO-1
    z2 = dcast(mapping_table[,.(UNIPROT, GO, value = 1)], UNIPROT ~ GO, fill = 0, drop = F)[,UNIPROT := NULL][,terms_to_compare, with=F]
    # sanity check: print(str(z2))
    # calculate what are the unique comparisons of columns
    combinations = t(combs(colnames(z2),2))
    # sanity check: print(as.data.table(combinations))
    if(parallel == T){
        # create cluster
        cl <- makeCluster(detectCores()-1)  
        # get library support needed to run the code
        clusterEvalQ(cl,library(data.table))
        # put objects in place that might be needed for the code
        # clusterExport(cl,c("combinations", "z2"))
        # parSapply
        dist = t(parSapply(cl = cl,as.data.table(combinations), function(x) kappa_score(z2[,c(x[1],x[2]),with = F])))
        # attach labels to distance calculations
        dist = cbind(as.data.table(dist), as.data.table(t(combinations)))
        #stop the cluster
        stopCluster(cl)
    }
    if(parallel == F){
        dist = t(sapply(as.data.table(combinations), function(x) kappa_score(z2[,c(x[1],x[2]),with = F])))
        dist = cbind(as.data.table(dist), as.data.table(t(combinations)))
    }
    colnames(dist) = c("kappa_score", "kappa_error", "GO1", "GO2")
    
    # convert the table to similarity matrix
    # generate the table with reverse comparisons (redundant) to make converting to avoid NAs in the symmetric similarity matrix
    dist_temp = unique(rbind(dist,dist[,.(kappa_score,kappa_error, GO1 = GO2, GO2 = GO1)]))
    # transform to matrix to find rownames
    dist2 = as.matrix(dcast(dist_temp[,.(GO1,GO2, kappa_score)], GO1 ~ GO2))
    rownames_dist2 = dist2[,"GO1"]
    # transform to matrix to get resulting similarity matrix and name the row correctly
    dist2 = as.matrix(dcast(dist_temp[,.(GO1,GO2, kappa_score)], GO1 ~ GO2)[,GO1 := NULL])
    rownames(dist2) = rownames_dist2
    # sort by row names and by column names
    dist2 = dist2[sort(rownames(dist2)), sort(colnames(dist2))]
    diag(dist2) = 1
    return(list(similarity_matrix = dist2, kappa_score_table = dist, kappa_score_table_redundant = dist_temp))
}