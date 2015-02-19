# black.miRNAmRNA.NSCLC.networkEnrichment

This `R` package contains functions and a vignette describing the generation of a STRING based protein-protein interaction (PPI) network for miRNA and mRNA data for NSCLC cells for the publication:

W Wu, RM Flight, MJ Krentz, B Kulengowski, H-F Li, HNB Moseley, A Stromberg and EP Black, *Interacting miRNA and mRNA genes identify potential therapeutic targets for NSCLC*, sumbitted Feb 2015

## Installation

If you want to be able to install this package directly and use it without modification to the *vignette*, then you are best doing a clone and then install using `devtools`, as the STRING database files are installed in `inst/extdata/STRING` so that there is no confusion about their location.

```
git clone https://github.com/rmflight/black_mirna_mrna_nsclc_network_enrichment.git
cd black_mirna_mrna_nsclc_network_enrichment/inst/extdata/STRING
wget http://string-db.org/newstring_download/protein.links.detailed.v9.1/9606.protein.links.detailed.v9.1.txt.gz
wget http://string.uzh.ch/permanent/string/9_1/protein_aliases/9606__protein_aliases_tf.tsv.gz

cd ../../../..
R CMD INSTALL black_mirna_mrna_nsclc_network_enrichment
```

To run the enrichment vignette, you will also need to install the following packages:

in R

```
library(BiocInstaller)
biocLite(c("graph", "org.Hs.eg.db", "hgu133plus2.db", "GO.db"))
```

If you want to re-generate the network figure, you will also need to install `Cytoscape 2.8` with the `CytoscapeRPC` plugin, and the `RCytoscape` Bioconductor package:

```
biocLite("RCytoscape")
```

## Original Output

All of the original tables and figures generated for the manuscript are in the `inst/extdata/output` directory

## Re-running Vignettes

If you want to re-run the vignettes for the network generation and enrichment, you can generate them all by using:

```
library(devtools)
build_vignettes()
```

Alternatively, you can use `knitr` directly:

```
library(knitr)
knit("vignettes/networkAnalysis.Rmd")
knit("vignettes/figureGeneration.Rmd")
```

## License

```
Licensed under the MIT License. This software is provided as-is, with no implied warranty.
```



