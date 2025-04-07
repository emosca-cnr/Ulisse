<img src="vignettes/images/logo.png" width="250">

# Ulisse - An R package to go beyond the boundaries of gene-sets

The understanding of how gene-related molecular alterations translate into pathological phenotypes is a major challenge in life sciences. 
Here, we address the challenge of assessing the possible alteration, with respect to a reference condition, of intra- and inter-cellular molecular interactions among sets of genes, which are intended to represent intra-cellular or cellular phenotypes. 
We provide a means to screen the alteration of intra-cellular pathway crosstalks and derive a map of the altered communications among pathways, which complements pathway enrichment analysis. 
Ulisse can also be used to reconstruct a cell-cell communication network between cell types/clusters. These two analyses (intra- and inter-cellular) can be combined to obtain integrated pathways of interactions that associate cell-cell communications with intracellular states. 
Further, we provide a score and a statistical assessment of the altered interactions, based on multiple empirical null models for networks. Lastly, we extract the key genes that take part in the altered interactions. 

Ulisse provide the tools to perform:

- Cross-talk analysis: intra-cellular and inter-cellular
- gene classification analysis: to reconstruct the role of the genes in the cross-talk network obtained
- Integrated cross-talk analysis: to integrate the communications between specific cells types/clusters to cellular mechanisms.

Typical application of Ulisse includes:

- Intra-cellular cross-talk analysis of omics data obtained from bulk or single-cell samples
- Inter-cellular analysis between clusters or cell-type in single-cell samples

Source code: https://github.com/emosca-cnr/Ulisse

# Installation

Ulisse requires R >= 4.0.0, and some GitHub and Bioconductor packages.

To successfully install Ulisse firstly run 

```{r, include=TRUE, eval=FALSE}
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(c("BiocParallel", "ComplexHeatmap", "DOSE", "KEGGREST", "qvalue"))

devtools::install_github("emosca-cnr/NPATools", build_vignettes = TRUE)
```

The other dependencies, if missing, should be automatically installed using the following command:

```{r, include=TRUE, eval=FALSE}
devtools::install_github("emosca-cnr/Ulisse", build_vignettes = TRUE)
```

Contacts:

- [Ettore Mosca](https://www.itb.cnr.it/en/institute/staff/ettore-mosca), Bioinformatics Lab, CNR-ITB
