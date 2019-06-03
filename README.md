# GSEA-Fisher: Gene Set Enrichment Analysis
Performs GeneSet Enrichment Analysis (GSEA) based on one-tail Fisher's Exact Test (Hypergeometric test). Implemented in R.
<br/><br/>

## *get.enrichment.test*

#### Decription:

Performs GeneSet Enrichment Analysis (GSEA) based on one-tail Fisher's Exact Test

#### Usage
 
```sh
get.enrichment.test(file.query, batch, file.background, dir.gmt, threshold, dir.out)
```

#### Arguments

- `file.query` : Path of the file containing the list of genes for which enrichment test is to be performed. One Gene per line
- `batch` : Name of the output folder. This folder will be under the folder `dir.out`
- `file.background` : Path of the file containing the list of background genes. One Gene per line
- `dir.gmt` : Path of the folder containing the genesets (.gmt) files
- `threshold` : Threshold of FDR p-value cutoff. The genesets with FDR less that **threshold** will be selected
- `dir.out` : Path of the output folder. Defaults to ''enrichment'' folder


#### Value 

The output file is a tab separated table containing following columns:

- `Category` :  Name of the Geneset
- `pvalue` : P-value of fisher's exact test
- `fdr` : False Discovery Rate (FDR) adjusted p-value, uses ''Benjamini Hochberg'' method for p-value correction
- `overlap.percent` : Percentage of query genes that are also found in the enriched geneset
- `overlap.genes` : List of query genes that are also found in the enriched geneset
- `Description` : Description of the geneset
<br/><br/>



## *get.fisher.exact.test*

#### Decription:

Performs one-tail Fisher's Exact Test (Hypergeometric test)

#### Usage
 
```sh
get.fisher.exact.test(dat.genesets, genes.queryset, genes.refset, ct)
```


#### Arguments

- `dat.genesets` : Data frame containing genesets (output file generated using **parseGMT()** function)
- `genes.queryset` : List of genes for which enrichment test is to be performed
- `genes.refset` : List of background set of genes
- `ct` : Threshold of FDR p-value cutoff. The genesets with FDR less that **ct** will be selected
<br/><br/>




## *parseGMT*

#### Decription:

Parse GMT file in appropriate format for enrichment test

#### Usage
 
```sh
parseGMT(dir.gmt)
```

#### Arguments

- `dir.gmt` : path of folder under which genesets (.gmt) files are stored


#### Value 

The output file is a tab separated table containing following columns:

- `Category` :  Name of the Geneset
- `Genesets` : List of genes in the geneset each seperated by a column (':')
- `Description` : Description of the geneset
<br/><br/>


## *How to run GSEA-Fisher*

Load the R script `run_GSEAfisher.R`

```{r echo=TRUE, message=FALSE, warning=FALSE}
source("run_GSEAfisher.R")
```

Parse GMT files. To be used only once for a file. If the GMT files are already parsed, skip this step.

```{r}
parseGMT(gmt.name="genesets/Msigdb")
```

Now run the enrichment test
```{r}
# Defile Paths
dir.db <- file.path("genesets/Msigdb")
file.bg <- file.path(dir.wrk, "data/background_genelist_test.txt")
file.genelist <- file.path(dir.wrk, "data/genelist_test.txt")

# Perform Enrichment Test
get.enrichment.test(file.query=file.genelist, batch="test", file.background=file.bg, dir.gmt=dir.db, threshold=0.01)
```

The output files can be found under `enrichment/TEST` folder.


