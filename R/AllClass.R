#' A class to store HiC data and all related objects for differential interaction analysis.
#'
#' @slot contactMatrixList A list stores HiC contact matrices
#' @slot d_height data frame of 3 columns storing xy coords and height/difference values
#' @slot persDiag data frame storing persistence function
#' @slot pLevelGrid a list whose elements contain different persistence levels of
#' \code{testingTree}
#' @slot path path to the working folder
#' @slot testingTree a list whose elements are nodes of \code{testingTree}
#' @slot checkTree a logical list to perform hierarchical testing
#' @slot hic_diff_result a matrix containing location of cell being tested with its p-values
#' @slot excluded  a matrix containing locations of cell not being tested
#'
#' @return
#' @export
#'
#' @examples
setClass("treeHiCDataSet",
         representation = representation(contactMatrixList = "list", d_height = "data.frame",
                                         persDiag = "data.frame", pLevelGrid = "list",
                                         path = "character", testingTree = "list",
                                         checkTree = "list", hic_diff_result = "matrix",
                                         excluded = "data.frame"),
         prototype = prototype(contactMatrixList = list())
)
