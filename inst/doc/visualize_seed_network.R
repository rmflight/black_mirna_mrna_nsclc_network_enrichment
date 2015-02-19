## ----load_packages, message=FALSE, warning=FALSE-------------------------
library(graph)
library(org.Hs.eg.db)
library(RCytoscape)

## ----load_data-----------------------------------------------------------
data(seedMRNA)
seedMRNA

data(PPI_seed)
head(PPI_seed)
nrow(PPI_seed)

## ----create_graph, warning=FALSE-----------------------------------------
seed_graph <- PPI_to_graph(PPI_seed, org.Hs.eg.db)

seedInfo <- seed_graph$info
seedGraph <- seed_graph$graph
seedGraph
head(seedInfo)

## ----collapse_complexes--------------------------------------------------
collapseString <- c("CNOT", "EIF3", "EXOC", "RPL", "RPS", "PPP1R", "PPP2R", "PPP3", "PPP4")
collapse_graph <- collapseNodes(seedGraph, seedInfo, seedMRNA$SYMBOL, collapseString)
collapseGraph <- collapse_graph$graph
collapseInfo <- collapse_graph$info

## ----setupRunningCytoscape, echo=FALSE, error=FALSE, message=FALSE-------
runCy <- FALSE

try({
  tmp <- RCytoscape::CytoscapeWindow("tmp")
  runCy <- TRUE
  invisible(RCytoscape::deleteWindow(tmp))
})

if (!runCy){
  print("No connection to Cytoscape available, subsequent visualizations were not run")
}

## ----visGraph, eval=runCy------------------------------------------------
cw <- CytoscapeWindow("collapsed", collapseGraph)
displayGraph(cw)
seedNodes <- collapseInfo$nodeID[collapseInfo$SYMBOL %in% seedMRNA$SYMBOL]
otherNodes <- grep("^9606", collapseInfo$nodeID, value=TRUE)
otherNodes <- otherNodes[!(otherNodes %in% seedNodes)]

setLayoutProperties(cw, layout.name="kamada-kawai-noweight", list(anticollisionStrength=4000, rest_length=100))

layoutNetwork(cw, "kamada-kawai-noweight")
redraw(cw)
setNodeLabelDirect(cw, nodes(collapseGraph), nodeData(collapseGraph, nodes(collapseGraph), "DISPLAY"))
redraw(cw)
setNodeColorDirect(cw, seedNodes, "#F5C46F")
redraw(cw)
setNodeColorDirect(cw, otherNodes, "#6FA0F5")
redraw(cw)
setDefaultBackgroundColor(cw, "#FFFFFF")
redraw(cw)
setDefaultNodeSize(cw, 60)
redraw(cw)
setDefaultEdgeColor(cw, "#CCCCCC")
redraw(cw)
setDefaultNodeFontSize(cw, "20")
redraw(cw)

