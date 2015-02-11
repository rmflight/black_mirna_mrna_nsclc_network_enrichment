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
#' @example
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