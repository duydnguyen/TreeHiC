.evalDiffMat <- function(object, useLog2, include.zeros, return_excluded_mat) {
    msMatList <- lapply(object@contactMatrixList, function(mat) {TreeHiC::convertForVtk(hicMat = mat)} )
    nreps <- length(msMatList)/2
    matC1 <- matC2 <- matRes <- msMatList[[1]]
    excluded <- matrix()
    meanf1 <- meanf2 <- rep(0, dim(matC1)[1])
    ## if more than 1 replicates for each cond., taking mean
    for (i in 1:nreps) {
        meanf1 <- meanf1 + msMatList[[i]][, 'f']
        meanf2 <- meanf2 + msMatList[[nreps + i]][, 'f']
    }
    meanf1 <- meanf1/nreps
    meanf2 <- meanf2/nreps
    matC1[,'f'] <- meanf1
    matC2[,'f'] <- meanf2
    ## zeros <- sapply(msMatList, FUN = function(x) which(x[, 'f'] == 0))
    ## zeros <- unique(unlist(zeros))
    ## exclude cells with zero means in BOTH conditions
    temp_msMatList <- list('matC1' = matC1, 'matC2' = matC2)
    zeros <- sapply(temp_msMatList, FUN = function(x) which(x[, 'f'] == 0))
    zeros <- unique(intersect(zeros[[1]], zeros[[2]]))
    if (useLog2) {
        matRes[, 'f'] <- log2(matC2[, 'f'] + 1) - log2(matC1[, 'f'] + 1)
    }
    else {
        matRes[, 'f'] <- matC2[, 'f'] - matC1[, 'f']
    }
    # scale to [0, 1]
    matRes[, 'f'] <- (matRes[, 'f'] - min(matRes[, 'f']))/(max(matRes[, 'f']) - min(matRes[, 'f']))
    temp.matRes <- matRes
    ## unSelected  <- union(zeros, which(matRes[,'f'] == 0))
    unSelected  <- zeros
    if (!include.zeros) {
        ## exlude zeros from each contact maps, and zero from any pixel in matRes
        ## unSelected  <- union(zeros, which(matRes[,'f'] == 0))
        matRes <- matRes[-unSelected, ]
    }
    ## if (return_excluded_mat) {
    ##     matRes <- list('matRes' = matRes, 'excluded' = temp.matRes[unSelected, ])
    ## }
    object@d_height <- matRes
    object@excluded <- temp.matRes[unSelected, ]
    object
}

setMethod("evalDiffMat", signature("treeHiCDataSet"), .evalDiffMat)
