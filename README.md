## TFnetwork
`TFnetwork` is an R package for Differential-expression analysis (DE-analysis), construction of core transcription-factor(TF) regulatory network (based on `NetAct`: https://github.com/lusystemsbio/NetAct) and core TF-target network plotting, all in one. What you need to do is simple: just input the expression matrix, and you will get the differentially-expression genes(DEGs), TFs with altered activity, a heatmap for TF activity and expression, TF-target network file and network plot. Unlike `NetAct`, `TFnetwork` only supports 2 groups, such as experiment group vs negative control.

`TFnetwork` is licensed under the GNU General Public License v3.0 https://github.com/PhiliJ/TFnetwork/blob/main/LICENSE

### Author

Li Lei <2020320243@stu.cqmu.edu.cn>

### Ph.D. supervisor

Yajun Xie <yjxie@cqmu.edu.cn>

Qin Zhou <zhouqin@cqmu.edu.cn>

### Installation:

You can install `TFnetwork` like so:

``` r
library(devtools)
install_github("PhiliJ/TFnetwork", dependencies = T)
``` 

### Step 1 Preparation of expreesion matrix for later analyses
#### 1.1 Set up input and comparisons

For example, for treated group and NC group:

First thing first, you will need to prepare the expression table of genes FOR EACH SAMPLE. 
The first row of the table should be ENSEMBL or ENTREZID, and the second row
should be the row count of each gene, then name these table like below.

Then, you will need to define variables for expression data files, sample names, group assignments, and a comparison to be made between groups for later analyses.

``` r
files <- c("nc1.txt", "nc2.txt", "nc3.txt",
           "treat1.txt", "treat2.txt", "treat3.txt")
samplenames <- c("nc1", "nc2", "nc3", "treat1", "treat2", "treat3")
group <- as.factor(c("nc", "nc", "nc", "treat", "treat", "treat"))
compList <- c("treat-nc")
```

#### 1.2 Set up expression matrix and phenodata
After that, use `Cons_exp` to construct expression matrix (x) for later analyses, then add phenodata in x.

``` r
x <- Cons_exp(files, group, samplenames, compList)

phenoData <- new("AnnotatedDataFrame", data = data.frame(celltype = group))
rownames(phenoData) <- colnames(x$counts)
```

#### 1.3 Preprocess expression counts
Use `pre` to pre-process the expression data. It filters out lowly expressed genes and retrieves associated gene symbols for raw count data. 
This function is based on `Preprocess_counts` function in `NetAct` package, but allows users to adjust `min.count`, `min.total.count`, `large.n` and `min.prop`.
You can use `?pre` for details.

``` r
counts <- pre(counts = x$counts, groups = group, mouse = TRUE)
```


### Step 2 DE-analysis for DEGs

You can use either `DESeq2` or `limma` for DE-analysis.

#### 2.1 DE-analysis using DESeq2

Use `DiffDESeq2` after setting cutoff values of logFC (logfc), P-value (pval) and using adjusted P or not (padj). You can use `?DiffDESeq2` for details.

``` r
DErslt = DiffDESeq2 (counts = counts, phenodata, complist = compList, 
                     logfc = 3, pval = 0.001, padj = TRUE)
                     
write.table(DErslt$table, file = "DESeq2_result.txt", sep = "\t")
```

#### 2.2 DE-analysis using limma

Use `Difflimma` after setting cutoff values of logFC (logfc), P-value (pval) and using adjusted P or not (padj). You can use `?Difflimma` for details.

``` r
DErslt = Difflimma (counts = counts, phenodata = phenoData, complist = compList, 
                     logfc = 3, pval = 0.001, padj = TRUE)
                     
write.table(DErslt$table, file = "limma_result.txt", sep = "\t", row.names = FALSE)
```


### Step 3 Identify significantly enriched TFs using NetAct
#### 3.1 TF GSEA

First, you need to use the command `data("mDB")` to call the `mDB` dataset for mouse
or `data("hDB")` to call the `hDB` dataset for humans in R.

``` r
data("mDB")
```

Then, use `TF_Selection_p` to do the GSEA and select the TF with P-value as cutoff.
If you want to use q-value as cutoff instead, just use `TF_Selection` in NetAct.
You can use `?TF_Selection_p` for details.

``` r
gsearslts <- TF_Selection_p(GSDB = mDB, DErslt = DErslt, minSize = 5, nperm = 1000,
                            pval = 0.05, compList = compList, method = "binary",
                            nameFile = "TF_GSEA")

write.table(gsearslts$GSEArslt, file = "TF_GSEA.txt", sep = "\t", row.names = FALSE)

tfs <- gsearslts$tfs
tfs
```

In case where more or fewer TFs are needed, you can use `Reselect_TFs_p` to re-select the TF with P-value as cutoff.
If you want to use q-value as cutoff instead, just use `Reselect_TFs` in NetAct.
You can use `?Reselect_TFs_p` for details.

``` r
tfs <- Reselect_TFs_p(GSEArslt = gsearslts$GSEArslt, pval = 0.01)
tfs
```




