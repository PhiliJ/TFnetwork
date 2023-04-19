#' @title TF selection with P-value
#' @description Identifying enriched TFs using Gene Set Enrichment Analysis (GSEA) with P-value as the cutoff. The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param GSDB gene set database (a list of gene sets, each of which is comprised of a vector genes)
#' @param DErslt DEG results
#' @param minSize the minimum number of overlapping genes required for each gene set (a gene set filtering parameter, default: 5)
#' @param nperm the number of gene list permutations (default: 1000)
#' @param method fast: fgsea; R: r implementation of GSEA with a new permutation method; binary: R/C++ implementation for fast speed
#' @param pval p-value cutoff (default: 0.05)
#' @param compList a vector of comparisons, it needs to be consistent with DErslt from MicroDegs, RNAseqDegs_limma, and RNAseqDegs_DESeq.
#'                  GSEA is applied to each comparison 
#' @param ntop the number of top genes (selection by the top genes) (default: NULL, no selection by the top genes)
#' @param nameFile file name to save the GSEA results (default: NULL, no output to a file). 
#'                  The saved results can be reused later to adjust the TF selection parameters
#' @return a list of results: 
#'         GSEArslt: a dataframe of GSEA results (see TF_GSEA).
#'         tfs: a vector of selected TFs.
#' @export
TF_Selection_p = function(GSDB, DErslt, minSize=5, nperm = 5000, method = "binary", pval = 0.05, 
                          compList = NULL, ntop = NULL, nameFile = NULL) {
  
  if(is.null(compList) | length(compList) == 1) {
    gsearslt <- TF_GSEA(GSDB = GSDB, DErslt = DErslt, minSize = minSize, nperm = nperm, method = method)
    tfs <- as.character(gsearslt$tf[gsearslt$pvals < pval[1]])
    if(!is.null(ntop)){
      if(length(tfs) > ntop){
        tfs = tfs[1:ntop]
      }
    }
    tfs = sort(tfs)
    output <- list(GSEArslt = list(gsearslt), tfs = tfs)
  }else{
    if ((length(pval) != 1) & (length(pval) != length(compList)))  {
      stop("pvalue length must be equal to 1 OR the number of comparisons")
    }
    if (length(pval) == 1) {
      pval <- rep(pval, length(compList))
    }
    gsea_results_list <- vector(mode = "list", length = length(compList))
    names(gsea_results_list) <- compList
    tfs = character()
    for (i in 1:length(compList)) {
      gsea_results_list[[i]] <- TF_GSEA(GSDB = GSDB, DErslt = DErslt[[i]], minSize = minSize, nperm = nperm, method = method)
      tfs_tmp <- as.character(gsea_results_list[[i]]$tf[gsea_results_list[[i]]$pvals < pval[i]])
      if(!is.null(ntop)){
        if(length(tfs_tmp) > ntop){
          tfs_tmp = tfs_tmp[1:ntop]
        }
      }
      tfs = c(tfs, tfs_tmp)
    }
    tfs = sort(unique(tfs))
    output <- list(GSEArslt = gsea_results_list, tfs = tfs)
  }
  
  if (!is.null(nameFile)) {
    saveRDS(output, file = paste0(nameFile, ".RDS" ))
  }
  return(output)
}

#' @title TF reselection with P-value
#' @description Reselecting TFs using gene set enrichement analysis (GSEA) using an adjusted set of parameters (work together with TF_Selection). The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param GSEArslt GSEA results from TF_Selection
#' @param pval p-value cutoff (default: 0.05)
#' @param combine_TFs whether combine selected TFs from multiple comparisons or not (default: TRUE)
#' @param ntop the number of top genes (selection by the top genes) (default: NULL, no selection by the top genes)
#' @return tfs: a vector of selected TFs
#' @export
Reselect_TFs_p = function(GSEArslt, pval = 0.05, combine_TFs = TRUE, ntop = NULL) {
  if (length(GSEArslt) == 1) {
    tfs <- as.character(GSEArslt[[1]]$tf[which(GSEArslt[[1]]$pvals < pval[1])])
    if(!is.null(ntop)){
      if(length(tfs) > ntop){
        tfs = tfs[1:ntop]
      }
    }
    tfs = sort(tfs)
  } else {
    if ((length(pval) != 1) & (length(pval) != length(GSEArslt)))  {
      stop("pvalue length must be equal to 1 OR the number of comparisons")
    }
    if (length(pval) == 1) {
      pval <- rep(pval, length(GSEArslt))
    }
    tfs = vector(mode = "list", length = length(GSEArslt))
    for (i in 1:length(GSEArslt)) {
      tfs_tmp <- as.character(GSEArslt[[i]]$tf[GSEArslt[[i]]$pvals < pval[i]])
      if(!is.null(ntop)){
        if(length(tfs_tmp) > ntop){
          tfs_tmp = tfs_tmp[1:ntop]
        }
      }
      tfs[[i]] = sort(tfs_tmp)
    }
    if (combine_TFs) {
      tfs = sort(unique(unlist(tfs)))
    } else {
      names(tfs) <- names(GSEArslt)
    }
  }
  return(tfs)
}

#' @title Row normalization (standardization)
#' @description Row normalization of gene expression matrix. The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param data gene expression matrix
#' @export
#' @return norm_data: standardized gene expression matrix
row_norm = function(data){
  row_mean = apply(data, 1, mean)
  row_sd = apply(data, 1, sd)
  norm_data = apply(data, 2, function(x) (x - row_mean)/(row_sd))
  return(norm_data)
}

#' @title TF reselection with P-value
#' @description Plotting heatmap of TF activity & TF gene expresion. The code is primarily based on the NetAct package, which can be found at https://github.com/lusystemsbio/NetAct.
#' @param new_activity a numeric matrix or data.frame containing the activity values of transcription factors (TFs)
#' @param eset an ExpressionSet or numeric matrix or data.frame containing the expression values of genes
#' @param activity_range a vector with three elements specifying the range of the activity values of TFs for color scaling
#' @param activity_color a vector with three colors specifying the color scheme for the activity heatmap
#' @param exp_range a vector with three elements specifying the range of the expression values of genes for color scaling
#' @param exp_color a vector with three colors specifying the color scheme for the expression heatmap
#' @return A heatmap object generated from the given TF activity and gene expression data
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap 
#' @export
TF_heatmap = function(new_activity, eset,
                      activity_range = c(-1, 0, 1),
                      activity_color = c("#4090BD", "#F5F5F5", "#E85355"),
                      exp_range = c(-2, 0, 2),
                      exp_color = c("#208421","#f0e9db50","#b40001")
){
  #    require(circlize)
  if (is(eset, "ExpressionSet")){
    data = exprs(eset)
  }else{data = eset}
  H1 = Heatmap(row_norm(new_activity), 
               col = colorRamp2(activity_range, activity_color), 
               rect_gp = gpar(col = "#F5F5F5", lwd = 0.5),
               cluster_columns = F, cluster_rows = T, 
               show_row_dend = FALSE, name = "Activity", column_title = "TF Activity")
  gs = rownames(new_activity)
  gc = colnames(new_activity)
  H2 = Heatmap(row_norm(data[gs, gc]), 
               col = colorRamp2(exp_range, exp_color),
               rect_gp = gpar(col = "#f0e9db", lwd = 0.5),
               cluster_columns = F, cluster_rows = T, 
               show_row_dend = FALSE,name = "Expression", column_title = "TF Expression")
  H1 + H2
  
}
