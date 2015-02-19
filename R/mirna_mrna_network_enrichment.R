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

#' get PPI network from STRING
#' 
#' Given a set of gene identifiers and a STRING data.frame, get the interactors within
#' a single edge of the \emph{seed} genes, and the set of all edges between the genes
#' 
#' @param seedGenes the genes to seed the PPI
#' @param stringDB the STRING data.frame to use
#' 
#' @export
#' @return data.frame
getPPI <- function(seedGenes, stringDB){
  hasSeed <- which((stringData$protein1 %in% seedGenes) | (stringData$protein2 %in% seedGenes))
  useNodes <- unique(c(stringDB$protein1[hasSeed], stringDB$protein2[hasSeed]))
  
  hasP1 <- stringDB$protein1 %in% useNodes
  hasP2 <- stringDB$protein2 %in% useNodes
  
  useEdges <- hasP1 & hasP2
  
  stringEdges <- stringDB[useEdges,]
  
  return(stringEdges)
  
}


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

#' @name PPI_seed.RData
#' @title PPI_seed.RData
#' @docType data
#' @source filtered STRING PPI network with edges > 400, seeded from EIF4B, RAB8A, and PPPC1B
NULL

#' @name seedMRNA.RData
#' @title seedMRNA.RData
#' @docType data
#' @source the \code{outMRNA_mapped} object for the seed genes of EIF4B, RAB8A, and PPPC1B
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


#' graph from PPI
#' 
#' generates a graph from the PPI network suitable for visualization
#' 
#' @param PPI the STRING interaction network
#' @param orgDB an organism database to add data to the network
#' @export 
#' @return list
#' @import graph
PPI_to_graph <- function(PPI, orgDB){
  ppiNodes <- unique(c(PPI$protein1, PPI$protein2))
  ppiGraph <- graphNEL(nodes=ppiNodes, edgemode="directed")
  ppiGraph <- addEdge(PPI$protein1, PPI$protein2, ppiGraph)

  trimNodes <- substring(ppiNodes, 6, 20)
  ppiInfo <- select(orgDB, trimNodes, c("SYMBOL", "GENENAME", "ENTREZID"), keytype="ENSEMBLPROT")
  ppiInfo$nodeID <- paste("9606.", ppiInfo$ENSEMBLPROT, sep="")

  for (iName in c("ENSEMBLPROT", "SYMBOL", "GENENAME", "ENTREZID")){
    nodeDataDefaults(ppiGraph, iName) <- ""
    attr(nodeDataDefaults(ppiGraph, iName), "class") <- "STRING"
    nodeData(ppiGraph, ppiInfo$nodeID, iName) <- ppiInfo[,iName]
  }
  
  return(list(info = ppiInfo, graph = ppiGraph))
}


#' collapse nodes and info
#' 
#' sometimes our graph has a lot of rather duplicated nodes that we want to collapse for visualization.
#' thats what this does. Given a graph, the info about the nodes, and the strings to use to find things to
#' collapse together, returns a new graph and info entries.
#' 
#' @param ppiGraph the graph of interactions
#' @param ppiInfo the information about each node
#' @param ppiSeed the nodes that were originally used to seed the PPI graph
#' @param collapseString the strings used to find those things to collapse
#' 
#' @export
#' @return list
collapseNodes <- function(ppiGraph, ppiInfo, ppiSeed, collapseString){
  dupEntry <- duplicated(ppiInfo$nodeID)
  ppiInfo <- ppiInfo[!dupEntry,]
  
  ppiInfo$collapseID <- ppiInfo$SYMBOL
  ppiInfo <- ppiInfo[!(is.na(ppiInfo$SYMBOL)),]
  
  for (checkCollapse in collapseString){
    whichCollapse <- grep(checkCollapse, ppiInfo$SYMBOL)
    ppiInfo$collapseID[whichCollapse] <- checkCollapse
  }
  
  ## ----collapseGraph Nodes-------------------------------------------------------
  edgemode(ppiGraph) <- "undirected"
  
  ppiMatrix <- as(ppiGraph, "matrix")
  
  ppiNew <- ppiMatrix
  ppiNew <- rbind(ppiNew, matrix(0, nrow=length(collapseSearch), ncol=ncol(ppiNew)))
  ppiNew <- cbind(ppiNew, matrix(0, nrow=nrow(ppiNew), ncol=length(collapseSearch)))
  rownames(ppiNew) <- c(rownames(ppiMatrix), collapseString)
  colnames(ppiNew) <- c(colnames(ppiMatrix), collapseString)
  
  for (iCollapse in collapseString){
    useID <- ppiInfo$nodeID[(ppiInfo$collapseID %in% iCollapse)]
    
    for (inID in useID){
      hasRow <- rownames(ppiNew)[ppiNew[,inID] == 1]
      ppiNew[hasRow, iCollapse] <- 1
      
      hasCol <- colnames(ppiNew)[ppiNew[inID,] == 1]
      ppiNew[iCollapse, hasCol] <- 1
    }
  }
  
  for (iCollapse in collapseString){
    ppiNew[iCollapse, iCollapse] <- 0
  }
  
  keepNodes <- c(ppiInfo$nodeID[!(ppiInfo$collapseID %in% collapseString)], collapseString)
  ppiNew <- ppiNew[keepNodes, keepNodes]
  
  ppiNewGraph <- as(ppiNew, "graphNEL")
  
  ## ----ppirabData----------------------------------------------------------
  ppiNewInfo <- ppiInfo[(ppiInfo$nodeID %in% nodes(ppiNewGraph)),]
  collapseData <- data.frame(ENSEMBLPROT=collapseString, SYMBOL=collapseString, GENENAME=collapseString, ENTREZID=collapseString, nodeID=collapseString, collapseID=collapseString, stringsAsFactors=FALSE)
  
  ppiNewInfo <- rbind(ppiNewInfo, collapseData)
  
  displayNodeLoc <- ppiNewInfo$SYMBOL %in% ppiSeed
  ppiNewInfo$DISPLAY <- ""
  ppiNewInfo[displayNodeLoc, "DISPLAY"] <- ppiNewInfo[displayNodeLoc, "SYMBOL"]
  
  for (iName in c("ENSEMBLPROT", "SYMBOL", "GENENAME", "ENTREZID", "DISPLAY")){
    nodeDataDefaults(ppiNewGraph, iName) <- ""
    attr(nodeDataDefaults(ppiNewGraph, iName), "class") <- "STRING"
    nodeData(ppiNewGraph, ppiNewInfo$nodeID, iName) <- ppiNewInfo[,iName]
  }
  ppiNewGraph <- initEdgeAttribute(ppiNewGraph, "weight", "numeric", 1)
  
  return(list(info = ppiNewInfo, oldInfo = ppiInfo, graph = ppiNewGraph))
}