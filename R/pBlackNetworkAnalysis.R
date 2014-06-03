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