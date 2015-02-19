## ----setup, message=FALSE, warning=FALSE---------------------------------
library(black.miRNAmRNA.NSCLC.networkEnrichment)
library(graph)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(GO.db)
library(xtable)

packageData <- system.file("extdata", package = "black.miRNAmRNA.NSCLC.networkEnrichment")
stringFileLoc <- file.path(packageData, "STRING")
inputDataLoc <- file.path(packageData, "input")

currLoc <- basename(getwd())
outputLoc <- switch(currLoc,
                    mirna_mrna_nsclc_network_enrichment = "inst/extdata/output",
                    vignettes = "../inst/extdata/output")

## ----loadData------------------------------------------------------------
mirData <- read.table(file.path(inputDataLoc, "compiledinteractiondata.csv"), header=TRUE, stringsAsFactors=FALSE, sep=",")
newmRNA <- sapply(seq(1, nrow(mirData)), function(x){
  useID <- mirData$mRNA[x]
  nID <- nchar(useID)
  substring(useID, 2, nID)
})
mirData$mRNA <- newmRNA

## ----linkDistribution----------------------------------------------------
splitMRNA <- split(mirData$miRNA, mirData$mRNA)
splitCount <- sapply(splitMRNA, length)
hist(splitCount)

## ----filterMRNA----------------------------------------------------------
nMin <- 5
sum(splitCount >= nMin)

useMRNA <- data.frame(gene=names(splitCount)[splitCount >= nMin], stringsAsFactors=FALSE)

## ----convertID-----------------------------------------------------------
outMRNA <- select(hgu133plus2.db, useMRNA$gene, c("SYMBOL", "GENENAME", "ENTREZID"), "PROBEID")

## ----getSTRINGID---------------------------------------------------------
stringIDFile <- file.path(stringFileLoc, "9606__protein_aliases_tf.tsv.gz")
outMRNA_mapped <- map2STRING(outMRNA, "SYMBOL", stringIDFile)
outMRNA_mapped$ENSEMBLPROT <- substring(outMRNA_mapped$protein_id, 6, 20)

# genes with larger expression values
seed_genes <- unique(outMRNA_mapped$SYMBOL)
outMRNA_mapped <- outMRNA_mapped[(outMRNA_mapped$SYMBOL %in% seed_genes),]

## ----stringData----------------------------------------------------------
stringData <- read.table(file.path(stringFileLoc, "9606.protein.links.detailed.v9.1.txt.gz"), header=TRUE, sep=" ", stringsAsFactors=FALSE)
head(stringData)

## ----filterData----------------------------------------------------------
combScore <- 400
minScores <- c(experimental=400, database=400, coexpression=400)
passesScore <- rep(FALSE, nrow(stringData))
invisible(lapply(names(minScores), function(evidence){
  passesScore <<- passesScore | (stringData[, evidence] >= minScores[evidence])
}))
passesScore <- passesScore & (stringData$combined_score >= combScore)
head(stringData[passesScore,])

## ----trueFilter----------------------------------------------------------
stringData <- stringData[passesScore,]

## ----getPPI_22-----------------------------------------------------------
PPI_22 <- getPPI(outMRNA_mapped$protein_id, stringData)

n_edges <- nrow(PPI_22)
uniq_genes <- unique(c(PPI_22$protein1, PPI_22$protein2))

## ----ensemblp_to_entrez--------------------------------------------------
library(GO.db)
all_proteins <- substring(unique(c(stringData$protein1, stringData$protein2)), 6)

protein_entrez <- select(org.Hs.eg.db, all_proteins, "ENTREZID", "ENSEMBLPROT")
protein_entrez <- protein_entrez[!(is.na(protein_entrez$ENTREZID)),]

gene_universe <- unique(protein_entrez$ENTREZID)
ensembl_diff <- substring(unique(c(PPI_22$protein1, PPI_22$protein2)), 6)
gene_diff <- unique(protein_entrez$ENTREZID[(protein_entrez$ENSEMBLPROT %in% ensembl_diff)])

gene_to_go <- select(org.Hs.eg.db, gene_universe, c("GOALL", "ONTOLOGYALL"), "ENTREZID")
gene_to_go <- gene_to_go[(gene_to_go$ONTOLOGYALL %in% "BP"), ]

go_annotation <- split(gene_to_go$ENTREZID, gene_to_go$GOALL)
go_annotation <- lapply(go_annotation, unique)

go_description <- Term(names(go_annotation))

## ----go_enrichment-------------------------------------------------------
go_annotation <- new("annotation",
                     annotation_to_feature = go_annotation,
                     description = go_description)

go_hypergeom <- new("hypergeom_features",
                    significant = gene_diff,
                    universe = gene_universe,
                    annotation = go_annotation)

sig_go_enrichment <- hypergeometric_feature(go_hypergeom)

sig_go_table <- sig_go_enrichment@annotation@stats
sig_go_table$description <- go_description[rownames(sig_go_table)]
write.table(sig_go_table, file = file.path(outputLoc, "go_22.tab"), row.names = TRUE, col.names = TRUE, sep = "\t")

## ----kegg_enrichment-----------------------------------------------------
data(kegg_annotation)
kegg_annotation <- new("annotation",
                       annotation_to_feature = kegg_annotation$annotation,
                       description = kegg_annotation$description)
kegg_hypergeom <- new("hypergeom_features",
                      significant = gene_diff,
                      universe = gene_universe,
                      annotation = kegg_annotation)

sig_kegg_enrichment <- hypergeometric_feature(kegg_hypergeom)
sig_kegg_table <- sig_kegg_enrichment@annotation@stats
sig_kegg_table$description <- kegg_annotation@description[rownames(sig_kegg_table)]
write.table(sig_kegg_table, file = file.path(outputLoc, "kegg_22.tab"), row.names = TRUE, col.names = TRUE, sep = "\t")

## ----visNetwork, eval=FALSE----------------------------------------------
#  seedGenes <- c("EIF4B", "RAB8A", "PPP1CB")
#  
#  seedMRNA <- outMRNA_mapped[(outMRNA_mapped$SYMBOL %in% seedGenes),]
#  
#  PPI_seed <- getPPI(seedMRNA$protein_id, stringData)
#  save(PPI_seed, file = "inst/data/PPI_seed.RData")
#  save(seedMRNA, file = "inst/data/seedMRNA.RData")

