## TFnetwork
`TFnetwork` is an R package for Differential-expression analysis (DE-analysis), construction of core transcription-factor(TF) regulatory network (based on `NetAct`: https://github.com/lusystemsbio/NetAct) and core TF-target network plotting, all in one. 

What you need to do is simple: just input the expression matrix, and you will get the differentially-expression genes(DEGs), TFs with altered activity, a heatmap for TF activity and expression, TF-target network file and network plot. 

Unlike `NetAct`, `TFnetwork` only supports 2 groups, such as experiment group vs negative control.

`TFnetwork` is licensed under the GNU General Public License v3.0 https://github.com/PhiliJ/TFnetwork/blob/main/LICENSE

### Author

Li Lei <2020320243@stu.cqmu.edu.cn>

Guozhi Zhao <2020320060@stu.cqmu.edu.cn>

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
If you want to use q-value as cutoff instead, just use `TF_Selection` in `NetAct`.
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
If you want to use q-value as cutoff instead, just use `Reselect_TFs` in `NetAct`.
You can use `?Reselect_TFs_p` for details.

``` r
tfs <- Reselect_TFs_p(GSEArslt = gsearslts$GSEArslt, pval = 0.01)
tfs
```

#### 3.2 Calculate the activities of the selected TFs
TF activities are calculated using `TF_Activity` in `NetAct`.

``` r
neweset <- Biobase::ExpressionSet(assayData = as.matrix(DErslt$e), phenoData = phenoData)
act.me <- TF_Activity(tfs, mDB, neweset, DErslt)
acts_mat = act.me$all_activities
```

Then, use `TF_heatmap` to draw the heatmap with TF activities and expressions.

``` r
pdf("TF_heatmap.pdf", width = 4, height = 9)
heatmap2 <- TF_heatmap(acts_mat, neweset)
heatmap2
dev.off()
```

![TF_heatmap](https://user-images.githubusercontent.com/39685949/233879806-08290313-1676-4e7a-b4bd-a06da0a45591.png)

### Step 4 TF-Target Network
#### 4.1 Construct TF-target network
To avoid having too many TFs in the network, you might consider selecting only the most relevant TFs to include in the network visualization.

``` r
tfs <- Reselect_TFs_p(GSEArslt = gsearslts$GSEArslt, pval = 0.00001)
tfs
```

Then, use `Cons_net` to construct the network. This will also return an '.txt' file contains TF-target, and you can use this to generate network plot
using other software like `Cytoscape`.

``` r
Cons_net(mouse = TRUE)
```

#### 4.2 Plot TF-target network
Use `plot_tfnetwork`. Red color stands for "up-regulated" while blue stands for "down-regulated". You may need to check the details using `?plot_tfnetwork`. 

``` r
plot_tfnetwork(net,
               pdf.name = "TFnet.pdf",
               nodeshape = "circle",
               pdf.height = 15,
               pdf.width = 15,
               edge.arrow.size = 0.3,
               vertex.label.family = "sans")
``` 

![TFnet](https://user-images.githubusercontent.com/39685949/233879380-716ef378-ff9b-46e4-890d-935738614e69.png)

Use `nodeshape = "sphere"` to draw a 3D network.

``` r
plot_tfnetwork(net,
               pdf.name = "TFnet.pdf",
               nodeshape = "sphere",
               pdf.height = 15,
               pdf.width = 15,
               edge.arrow.size = 0.3,
               vertex.label.family = "sans")
``` 
![TFnet_3d](https://user-images.githubusercontent.com/39685949/233880282-4bb26a40-011a-48b0-ad7d-447f59999054.png)


#### 4.3 Highlight nodes of interest
Use `plot_highlighted_tfnetwork`. `TFnetwork` will automatically determine whether your node of interest is a TF or a target. You may need to check the details using `?plot_highlighted_tfnetwork`.

For example, for TF like Gli2:

``` r
plot_highlighted_tfnetwork(network, highlight_node = "Gli2")
``` 
![TFnet_highlight_Gli2](https://user-images.githubusercontent.com/39685949/233879439-a2191bd4-07c3-4f26-bf41-318d841a1b5d.png)

And for target like Pax8:

``` r
plot_highlighted_tfnetwork(network, highlight_node = "Pax8")
``` 
![TFnet_highlight_Pax8](https://user-images.githubusercontent.com/39685949/233879481-6b7b03d8-c658-4839-9e7e-39f26307a9f3.png)


#### 4.4 Plot interactive network
Use `plot_interactive_tfnetwork`. This allows you to move the node and adjust the layout.

``` r
plot_interactive_tfnetwork(network,
                           canvas.width = 500, 
                           canvas.height = 500)
``` 
![TF interactive network](https://user-images.githubusercontent.com/39685949/233880041-2f1fd2fa-0efa-46da-8185-e56be2f27a72.png)

#### 4.4 Plot 3D network
Just for fun :)
``` r
plot_3d_tfnetwork(network)
``` 
![3d](https://user-images.githubusercontent.com/39685949/233880193-456b198d-25f4-46ab-910e-ed38f7ff2f53.png)
