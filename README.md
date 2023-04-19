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

### Step 1
#### Prepare the expreesion matrix for later analyses

For example, for treated group and NC group:

First thing first, you will need to prepare the expression table of genes. The first row of the table should be ENSEMBL or ENTREZID, and the second row
should be the row count of each gene, then name these table like below.
Then, you will need to define variables for expression data files, sample names, group assignments, and a comparison to be made between groups for later analyses.

``` r
files <- c("nc1.txt", "nc2.txt", "nc3.txt",
           "treat1.txt", "treat2.txt", "treat3.txt")
samplenames <- c("nc1", "nc2", "nc3", "treat1", "treat2", "treat3")
group <- as.factor(c("nc", "nc", "nc", "treat", "treat", "treat"))
compList <- c("treat-nc")
```
