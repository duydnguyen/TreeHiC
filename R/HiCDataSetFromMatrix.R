.hicDataSetFromMatrix <- function(object, contactMatrixList, colData, path) {
    cond_lab <- unique(colData$Condition)
    ids.order <- c(which(colData$Condition == cond_lab[1]), which(colData$Condition == cond_lab[2]))
    sample.ids.order <- as.character(colData[ids.order,]$SampleID)
    #TODO: file_path inputs
    file_paths <- as.character(colData[ids.order,]$files)
    object@contactMatrixList <- contactMatrixList[ids.order]
    names(object@contactMatrixList) <- sample.ids.order
    object@path <- path
    return(object)
}

setMethod("HiCDataSetFromMatrix", signature("treeHiCDataSet"), .hicDataSetFromMatrix)
