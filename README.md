This is a repository for the coding of the `read10xMultiome` function.
This procedure will eventually be placed somewhere in Bioconductor.

`read10xMultiome` take as input a path to a 10x Genomics Multiome experiment and
returns a `SingleCellExperiment` object with the RNA (or ATAC) as main assay and 
the ATAC (or RNA) as `altExp` assay.
Main and `altExp` assays are decided with the reference parameter.
Additionally, it provides the possibility to load the peaks annotations with the
`annotation` parameter and to store them in the `rowRanges` of the ATAC assay. 
This last described process is a drafted algorithm, so it is very slow at the
moment.