.selectPLevelGrid <- function(object, maxDepth) {
    PD_CellData <- readr::read_csv(paste0(object@path,"temp/PD-CellData.csv"),
                            col_types = cols(`Cell Type` = col_skip(),
                            PairIdentifier = col_skip(), PairType = col_skip()))
    id_zeros <- which(PD_CellData$Persistence == 0)
    PD_CellData$Persistence[id_zeros] <- 10^{-14}
    object@persDiag <- PD_CellData
    ## Main
    dat <- data.table::data.table('nExtremas' = dim(PD_CellData)[1]:1, 'persistence' = sort(PD_CellData$Persistence))
    dat <- dat[persistence <= 1]
    n <- dim(dat)[1]
    #### create lag function
    x_diff <- log10(dat$persistence[2:n] - dat$persistence[1:(n-1)]) #better: clear separation between [0,pLevel1] and [pLevel1, 1]
    # @ smooth
    x_diff <- TreeHiC::ma(x_diff,order= 3)
    dat_diff <- data.table::data.table('persistence' = dat$persistence[2:n], 'pers_diff' = x_diff)
    #### pick pLevel1
    lmax <- 1
    dat_diff_sort <- dat_diff[order(pers_diff, decreasing = TRUE),]
    pLevel_diff <- as.numeric(dat_diff_sort[lmax,'pers_diff'])
    pLevel1 <- as.numeric(dat_diff_sort[lmax,'persistence'])
    nExtrema <- as.numeric(dat[persistence == pLevel1, 'nExtremas'])
    pLevelGrid <- c()
    if (maxDepth == 0) { ## only use pLevel1
        pLevelGrid <- pLevel1
    } else {
        #### pick pLevel > pLevel1
        dat_sort <- dat[persistence <= pLevel1, ]
        dat_sort <- dat_sort[order(persistence, decreasing = TRUE),]
        ## probs <- unique(c(seq(0.5, 0.99, 10^{-2}),c(0.99, 1)))
        probs <- unique(c(seq(0.5, 0.999, 10^{-2}),c(0.999, 1)))
        pLevelGrid <- quantile(dat_sort$persistence, probs = probs)
        maxDepth <- min(length(pLevelGrid), maxDepth)
        pLevelGrid <- pLevelGrid[(length(pLevelGrid) - maxDepth):(length(pLevelGrid) - 1 )]
        pLevelGrid <- c(pLevelGrid, pLevel1)
    }
    # return
    object@pLevelGrid <- list('pLevelGrid' = pLevelGrid, 'nExtremaPLevel1' = nExtrema,
                              'dat_diff' = dat_diff, 'pLevel_diff' = pLevel_diff, 'dat' = dat)
    object
}

setMethod("selectPLevelGrid", signature("treeHiCDataSet"), .selectPLevelGrid)
