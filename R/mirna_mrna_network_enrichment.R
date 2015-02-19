#' package black.miRNAmRNA.NSCLC.networkEnrichment
#' 
#' Package used to generate STRING based protein-protein interactions and perform
#' network enrichment using KEGG pathways and GO biological process for the publication:
#' W Wu, RM Flight, MJ Krentz, B Kulengowski, H-F Li, HNB Moseley, A Stromberg and EP Black,
#' \emph{Interacting miRNA and mRNA genes identify potential therapeutic targets for NSCLC},
#' sumbitted Feb 2015.
#' 
#' @name black.miRNAmRNA.NSCLC.networkEnrichment
NULL

#' spit out community structure
#' 
#' @param commObject the community object returned by an igraph community calculation
#' @param decreasing sort the communities with largest first (TRUE, default) or last (FALSE)
#' @export
#' @return list of communities
communityAsList <- function(commObject, decreasing=TRUE){
  listCommunity <- split(commObject$names, commObject$membership)
  commSizes <- sapply(listCommunity, length)
  commOrder <- order(commSizes, decreasing=decreasing)
  listCommunity <- listCommunity[commOrder]
  return(listCommunity)
}

#' @name allKGML.RData
#' @title allKGML.RData
#' @docType data
#' @source downloaded using KEGGREST from KEGG on June 4, 2014
#' @details contains set of KGML xml data that P Blacks genes were annotated to as of June 4, 2014
NULL

#' @name kegg_annotation.RData
#' @title kegg_annotation.RData
#' @docType data
#' @source downloaded using KEGGREST on Feb 11, 2015
#' @details list with annotation contains the gene - pathway information, and pathway descriptions
#' @examples
#' \dontrun{
#' library(KEGGREST)
#' hsa_pathways <- keggLink("hsa", "pathway")
#' names(hsa_pathways) <- substring(names(hsa_pathways), 6)
#' hsa_pathways <- substring(hsa_pathways, 5)
#' kegg_annotation <- list(annotation = split(hsa_pathways, names(hsa_pathways)))
#' kegg_annotation$annotation <- lapply(kegg_annotation$annotation, function(x){names(x) <- NULL; x})
#' 
#' kegg_paths <- names(kegg_annotation$annotation)
#' kegg_desc <- keggList("pathway", "hsa")
#' names(kegg_desc) <- substring(names(kegg_desc), 6)
#' kegg_annotation$description <- kegg_desc[names(kegg_annotation$annotation)]
#' }
NULL

#' map to STRING
#' 
#' @param query what you want to match on
#' @param stringAliasFile the file with STRING aliases
#' @export
#' @return data.frame of translated queries
map2STRING <- function(queryData, queryCol, stringAliasFile, removeUnmapped=TRUE){
  aliasTable <- read.table(stringAliasFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", fill=TRUE)
  foundAlias <- aliasTable$alias %in% queryData[, queryCol]
  aliasTable <- aliasTable[foundAlias, c("protein_id", "alias")]
  
  outData <- merge(queryData, aliasTable, by.x=queryCol, by.y="alias", all=TRUE)
  
  if (removeUnmapped){
    naP <- is.na(outData$protein_id)
    outData <- outData[!naP,]
  }
  return(outData)
}

#' annotation class
#' 
#' @slot annotation_to_feature list where each named entry is the features annotated to it
#' @slot description named character vector
#' @slot link named character vector
#' 
#' @export
setClass("annotation",
         slots = list(annotation_to_feature = "list",
                      description = "character",
                      link = "character",
                      stats = "data.frame"))

#' hypergeom feature class
#' 
#' class to hold features undergoing hypergeometric enrichment
#' 
#' @slot significant the significant features
#' @slot universe all of the features measured
#' @slot annotation annotation object
#' 
#' @export
setClass("hypergeom_features",
         slots = list(significant = "ANY",
                      universe = "ANY",
                      annotation = "annotation"))

#' do hypergeometric enrichment
#' 
#' @param hypergeom_features a hypergeometric_features object
#' @param direction which direction to do the enrichment (over or under)
#' @export
#' @return hypergeometric_Feature
#' 
hypergeometric_feature <- function(hypergeom_features, direction = "over"){
  
  # cleanup the features and annotations (should be in separate function)
  hypergeom_features@universe <- unique(hypergeom_features@universe)
  
  tmp_annot_feature <- hypergeom_features@annotation@annotation_to_feature
  annotation_universe <- unique(unlist(tmp_annot_feature))
  
  hypergeom_features@universe <- intersect(hypergeom_features@universe, annotation_universe)
  tmp_annot_feature <- lapply(tmp_annot_feature, intersect, hypergeom_features@universe)
  
  n_feature <- sapply(tmp_annot_feature, length)
  keep_annot <- n_feature > 0
  
  tmp_annot_feature <- tmp_annot_feature[keep_annot]
  
  hypergeom_features@significant <- intersect(hypergeom_features@significant, hypergeom_features@universe)
  hypergeom_features@annotation@annotation_to_feature <- tmp_annot_feature
  
  
  # this probably needs its own function eventually
  if (length(hypergeom_features@annotation@description) != 0){
    hypergeom_features@annotation@description <- hypergeom_features@annotation@description[names(tmp_annot_feature)]
  }
  
  if (length(hypergeom_features@annotation@link) != 0){
    hypergeom_features@annotation@link <- hypergeom_features@annotation@link[names(tmp_annot_feature)]
  }
  
  # now get the counts annotated
  num_white_drawn <- sapply(hypergeom_features@annotation@annotation_to_feature, function(x) sum(hypergeom_features@significant %in% x))
  
  if (length(num_white_drawn) == 0){
    num_white_drawn <- 0
  }
  
  num_white <- Biobase:::listLen(hypergeom_features@annotation@annotation_to_feature)
  
  if (length(num_white) == 0){
    num_white <- 0
  }
  
  num_black <- length(hypergeom_features@universe) - num_white
  num_drawn <- length(hypergeom_features@significant)
  
  hyper_stats <- hypergeometric_basic(num_white, num_black, num_drawn, num_white_drawn, direction)
  hyper_stats <- as.data.frame(hyper_stats)
  hyper_stats <- hyper_stats[(order(hyper_stats$p, decreasing = FALSE)),]
  
  hyper_stats$counts <- num_white_drawn[rownames(hyper_stats)]
  
  hypergeom_features@annotation@stats <- hyper_stats
  
  hypergeom_features
  
}

#' do hypergeometric test
#' 
#' does a hypergeometric enrichment test
#' 
#' @param num_white number of white balls in urn
#' @param num_black number of black balls in urn
#' @param num_drawn number of balls taken from urn
#' @param num_white_drawn number of white balls taken from urn
#' @param direction which direction is the test
#' 
#' @export
#' @return list
hypergeometric_basic <- function(num_white, num_black, num_drawn, num_white_drawn, direction = "over"){
  n_2_1 <- num_white - num_white_drawn
  n_1_2 <- num_drawn - num_white_drawn
  n_2_2 <- num_black - n_1_2
  
  odds_ratio <- (num_white_drawn * n_2_2) / (n_1_2 * n_2_1)
  
  expected <- ((num_white_drawn + n_1_2) * (num_white_drawn + n_2_1)) / (num_white_drawn + n_1_2 + n_2_1 + n_2_2)
  
  p_values <- switch(direction,
                     over =  phyper(num_white_drawn - 1L, num_white, num_black, num_drawn, lower.tail = FALSE),
                     under = phyper(num_white_drawn, num_white, num_black, num_drawn, lower.tail = TRUE)
  )
  
  list(p = p_values, odds = odds_ratio, expected = expected)
}
