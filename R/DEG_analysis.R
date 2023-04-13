#' @title Construct expression matrix
#' @description Construct expression matrix
#' @param files pre-defined inputs
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
#' @description Preprocess expression coutns
#' @param mouse use mouse genome or not
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
#' @param counts Processed gene expression count data
#' @param phenodata pData that provides batch & experimental conditions
#' @param complist a vector of multiple comparisons in the format of contrasts in limma (e.g. c("A-B", "A-C", "B-C"))
#' @param pval p-value cutoff for DEG analysis (default: 0.05)
#' @param padj Use adjusted p-value or not (default: TRUE)
#' @return DEresult: a list of DEG results, including those for each single comparison and those for the overall comparison.
#'         Each DEG result is in the format of A list containing:
#'         table: table of DEG results.
#'         rank_vector: a vector of t-statistics for every gene.
#'         degs: a vector of gene symbols for DEGs.
#'         e: expression data (CPM).
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
#' @param counts Processed gene expression count data
#' @param phenodata pData that provides batch & experimental conditions
#' @param complist a vector of multiple comparisons in the format of contrasts in limma (e.g. c("A-B", "A-C", "B-C"))
#' @param logfc (optional) log fold change constraints for DEGs
#' @param pval P-value cutoff for DEG analysis (default: 0.05)
#' @param padj Use adjusted p-value or not (default: TRUE)
#' @return DEresult: a list of DEG results, including those for each single comparison and those for the overall comparison.
#'         Each DEG result is in the format of A list containing:
#'         table: table of DEG results.
#'         rank_vector: a vector of t-statistics for every gene.
#'         degs: a vector of gene symbols for DEGs.
#'         e: expression data (CPM).
#'         e_batch: batch corrected expression.
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

