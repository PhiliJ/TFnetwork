#' @title Construct expression matrix
#' @description Construct expression matrix for later analysis. The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param files A vector containing the file names of the data.
#' @param group A factor vector indicating the group membership of each sample.
#' @param samplenames A vector containing the names of the samples.
#' @param compList A character vector specifying the comparison of interest. 
#' @return x: expression matrix for later analysis
#' @export
#' @import edgeR
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi

Cons_exp <- function(files, group, samplenames, compList) {
  x <- edgeR::readDGE(files, columns=c(1,2))
  x$samples$group <- as.factor(group)
  colnames(x) <- samplenames
  return(x)
}

#' @title Preprocess expression coutns
#' @description Preprocess expression coutns. The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param counts A matrix or data frame containing the raw count data.
#' @param groups A factor vector indicating the group membership of each sample.
#' @param mouse A logical indicating whether the input data is from mouse or human.
#' @param min.count The minimum number of counts for a gene to be considered expressed in a sample.
#' @param min.total.count The minimum total count across all samples for a gene to be considered expressed.
#' @param large.n The threshold number of samples in which a gene must be expressed to avoid filtering due to low expression.
#' @param min.prop The minimum proportion of samples in which a gene must be expressed to avoid filtering due to low expression.
#' @return x$counts: processed count matrix
#' @export
#' @import edgeR
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
pre <- function(counts, groups, mouse = FALSE,
                min.count = 1, min.total.count = 2, large.n = 1, min.prop = 0.1) {
  #  require(edgeR)
  #  Create DGEList and Filter
  x <- DGEList(counts, group = groups)
  
  #  Gene Annotations and Filter Duplicate Symbols
  geneid <- rownames(counts)
  
  if(grepl("ENS", geneid[1])) {
    keytype = "ENSEMBL"
  } else {
    keytype = "ENTREZID"
  }
  
  if (mouse == TRUE) {
    genes <- AnnotationDbi::select(org.Mm.eg.db, keys = geneid, columns = "SYMBOL", keytype = keytype)
  } else {
    genes <- AnnotationDbi::select(org.Hs.eg.db, keys = geneid, columns = "SYMBOL", keytype = keytype)
  }
  
  dup_id = which(duplicated(genes[keytype]))
  if(length(dup_id) != 0){
    genes = genes[-dup_id,]
  }
  
  filter_na <- which(is.na(genes$SYMBOL))
  filter_dup <- which(duplicated(genes$SYMBOL))
  filter <- unique(filter_na, filter_dup)
  if(length(filter) != 0){
    x$counts <- x$counts[-filter, ]
    rownames(x$counts) <- genes$SYMBOL[-filter]
  }
  
  #Filter Low Expression
  keep.exprs <- filterByExpr(x, group=groups, 
                             min.count = min.count, min.total.count = min.total.count, 
                             large.n = large.n, min.prop = min.prop)
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  write.table(counts, file = "Exp_counts.txt", sep = "\t")
  return(x$counts)
}

#' @title DEG Analysis of RNA-seq Data using DESeq with P value
#' @description This function performs differential expression analysis using the DESeq2 package.
#' @param counts A count matrix where rows are genes and columns are samples.
#' @param phenodata A data frame containing the phenotype data for each sample.
#' @param complist A character vector containing the names of the samples to compare.
#' @param logfc The log2 fold change cutoff for differential expression. Default is 0.
#' @param pval The adjusted p-value or p-value cutoff for differential expression. Default is 0.05.
#' @param padj Logical indicating whether to use adjusted p-values or not. Default is TRUE.
#' @import DESeq2
#' @return DEresult: A list containing the following items:
#' table: A data frame with the results of the differential expression analysis.
#' rank_vector: A vector containing the rank of each gene based on the effect size.
#' degs: A character vector containing the names of the differentially expressed genes.
#' e: A matrix of normalized expression values.
#' up: A character vector containing the names of the upregulated genes.
#' down: A character vector containing the names of the downregulated genes.
#' nosig: A character vector containing the names of the non-differentially expressed genes.
#' @export
#' @import DESeq2
DiffDESeq2 = function(counts, phenodata, complist,logfc = 0, pval = 0.05, padj = TRUE ) {
  
  phenodata = data.frame(celltype = group, row.names = samplenames)
  celltype = phenodata$celltype
  dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
  dds <- estimateSizeFactors(dds)
  e <- counts(dds, normalized = TRUE)
  
  dds = DESeq(dds)
  if (padj) {
    DEtable = results(dds, alpha = pval)
    DEtable = data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$padj), ]
    degs = rownames(DEtable[abs(DEtable$log2FoldChange) > logfc & DEtable$padj < pval, ])
  } else {
    DEtable = results(dds, alpha = pval)
    DEtable = data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    degs = rownames(DEtable[abs(DEtable$log2FoldChange) > logfc & DEtable$pvalue < pval, ])
  }
  
  # 计算高低差异基因数目
  upregulated = rownames(DEtable[DEtable$log2FoldChange > logfc & (if (padj) DEtable$padj < pval else DEtable$pvalue < pval), ])
  downregulated = rownames(DEtable[DEtable$log2FoldChange < -logfc & (if (padj) DEtable$padj < pval else DEtable$pvalue < pval), ])
  no_difference = rownames(DEtable[abs(DEtable$log2FoldChange) <= logfc | (if (padj) DEtable$padj >= pval else DEtable$pvalue >= pval), ])
  
  # 输出结果
  num_upregulated = length(upregulated)
  num_downregulated = length(downregulated)
  num_no_difference = length(no_difference)
  
  print(paste("Down:", num_downregulated))
  print(paste("NotSig:", num_no_difference))
  print(paste("Up:", num_upregulated))
  
  rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
  return(list(table = DEtable, rank_vector = rank_vector, degs = degs, e = e,
              up = upregulated, down = downregulated, nosig = no_difference))
}

#' @title DEG Analysis of RNA-seq Data using limma + Voom with P value
#' @description This function performs differential expression analysis using the limma package. The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param counts A matrix of raw counts.
#' @param phenodata A data frame containing phenotype information.
#' @param complist A list of comparisons to be performed.
#' @param logfc A numeric value indicating the log fold change threshold (default: 1).
#' @param pval A numeric value indicating the p-value threshold (default: 0.05).
#' @param padj A logical value indicating whether to adjust p-values (default: TRUE).
#' 
#' @return A list containing the following elements:
#' \item{table}{The differential expression table.}
#' \item{rank_vector}{A rank vector.}
#' \item{degs}{A list of differentially expressed genes.}
#' \item{e}{Normalized expression values.}
#' @export
#' @import edgeR
#' @import limma
Difflimma = function(counts, phenodata, complist, logfc = 1, pval = 0.05, padj = TRUE) {
  #  require(edgeR)
  #  require(limma)
  phenodata = data.frame(celltype = group, row.names = samplenames)
  celltype = phenodata$celltype
  levs = levels(celltype)
  dge = DGEList(counts=counts,group=celltype, genes = rownames(counts))
  dge <- calcNormFactors(dge)
  
  compList <- complist
  
  
  design = model.matrix(~ 0 + celltype)
  colnames(design) = levs
  
  contr.matrix = makeContrasts(contrasts = compList, levels = design)
  v <- voom(dge, design, plot = FALSE)
  vfit = lmFit(v, design)
  vfit = contrasts.fit(vfit, contrasts=contr.matrix)
  efit = eBayes(vfit)
  results = decideTests(efit, p.value = pval)
  
  #log fold change restriction
  if (!missing(logfc)) {
    
    efit <- treat(vfit, lfc  = logfc)
    if (padj) {
      results <- decideTests(efit, p.value = pval)
    } else {
      results <- decideTests(efit,adjust.method = "none", p.value = pval)
    } 
  }
  
  print(summary(results))
  
  data.rm.batch = NULL
  
  if (padj) {
    rslt = topTable(efit, coef=compList, number = Inf, sort.by = "P")
    rslt$padj = rslt$adj.P.Val
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< pval, ])
  } else {
    rslt = topTable(efit, coef=compList, number = Inf, sort.by = "P")
    rslt$padj = rslt$P.Value
    rank_vector = abs(rslt$t); names(rank_vector) = rownames(rslt)
    degs = rownames(rslt[rslt$padj< pval, ])
  } 
  
  return(list(table = rslt, rank_vector = rank_vector, degs = degs, 
              e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch)))
  
  
  
  DEtable = topTable(efit, coef = NULL, number = Inf, sort.by = "P")
  rank_vector = abs(DEtable$t); names(rank_vector) = rownames(DEtable)
  
  if (padj) {   
    DEtable$padj = DEtable$adj.P.Val
    degs = rownames(DEtable[DEtable$padj < pval, ])
  } else {
    DEtable$padj = DEtable$P.Value
    degs = rownames(DEtable[DEtable$padj < pval, ])
  }
  
  DEresult[["Overall"]] = list(table = DEtable, rank_vector = rank_vector, degs = degs,
                               e = as.data.frame(v$E), e_batch = as.data.frame(data.rm.batch))
  return(DEresult)
  
}

