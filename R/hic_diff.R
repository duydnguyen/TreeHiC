.hic_diff <- function(object, mat_pvals, use_adjusted_pvals, alpha,
                     batch_mode) {
    path_temp <- paste0(object@path,"temp/")
    pLevelGrid <- object@pLevelGrid[['pLevelGrid']]
    thres_names <- paste0('thresholdPLevel', 1:length(pLevelGrid))
    pathMSComplexPLevel1 <- paste0(path_temp, "MSComplexPLevel1.csv")
    path_vtk_coords <- paste0(path_temp, "vtk-coords.csv")
    pathThresholdPLevelList <- list()
    for (i in 1:length(pLevelGrid)) {
        pathThresholdPLevelList[[thres_names[i]]] <- paste0(path_temp, paste0(thres_names[i],'.csv'))
    }
    ## create testing tree
    testingTree <- TreeHiC::create_testing_tree(pathThresholdPLevelList = pathThresholdPLevelList,
                                               pathMSComplexPLevel1 = pathMSComplexPLevel1, d_height = object@d_height,
                                               batch_mode = batch_mode, path_vtk_coords = path_vtk_coords)
    testingTree <- TreeHiC::evalPvals(testingTree = testingTree, mat_pvals = mat_pvals)
    ## perform differential testing
    diffResultList <- TreeHiC::get_diff_hic(testingTree = testingTree, alpha = alpha, use_adjusted_pvals = use_adjusted_pvals)
    if (dim(diffResultList[['hic_diff_res']])[1] > 0) {
        object@hic_diff_result <- as.matrix(diffResultList[['hic_diff_res']])
    }
    object@checkTree <- diffResultList[['checkTree']]
    ## return object
    object@testingTree <- testingTree
    object
}

setMethod("hic_diff", signature("treeHiCDataSet"), .hic_diff)
