#' prepare HiC contact matrices for differential interaction. For this function, you should provide
#' count matrices, or txt files specified in the column information \code{colData}.
#'
#' @param object : a \code{treeHiCDataSet} object
#' @param contactMatrixList : a list of HiC contact matrices to perform differential interaction
#' @param colData : a table of sample information
#' @param path : path to the working folder
#'
#' @return a contact matrix list with respect to the column information \code{colData}
#' @export
#'
#' @examples
setGeneric("HiCDataSetFromMatrix", function(object, contactMatrixList = list(), colData,
                                            path) {
    standardGeneric("HiCDataSetFromMatrix")
})

#' create a normalized height function between two HiC matrices
#'
#' @param object : a \code{treeHiCDataSet} object
#' @param useLog2 : Logical, default to TRUE to use log2 of two matrices
#' @param include.zeros : should partial zero interactions be included? Defaul to TRUE
#'
#' @return a vector with range from 0 to 1 storing height difference
#' @export
#'
#' @examples
setGeneric("evalDiffMat", function(object, useLog2 = TRUE, include.zeros = TRUE, return_excluded_mat = FALSE) {
    standardGeneric("evalDiffMat")
})

#' pich the most persistent number of partitions and extremas
#'
#' @param object : a \code{treeHiCDataSet} object
#' @param maxDepth : maximum depth of testing tree
#' @param path : path to the working folder
#'
#' @return return the optimal persistent level
#' @export
#'
#' @examples
setGeneric("selectPLevelGrid", function(object, maxDepth = 0) {
    standardGeneric("selectPLevelGrid")
})


#' get cell ID (or vertexID vtk pipeline) for pLevel1 and pLevel2 with inputs
#' from Topological Toolkit (ttk)
#'
#' @param object : a \code{treeHiCDataSet} object
#' @param mat_pvals : a squared matrix of (raw) p-values to perform testing
#' @param use_adjusted_pvals : Logical, should adjusted p-values be used?
#' @param alpha : FDR control at \code{alpha} level
#'
#' @return add additional slots \code{testingTree, checkTree, hic_diff_result}
#' @export
#'
#' @examples
setGeneric("hic_diff", function(object, mat_pvals, use_adjusted_pvals = TRUE, alpha) {
    standardGeneric("hic_diff")
})
