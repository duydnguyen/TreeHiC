#' convert rep to for nx3 matrix: f(x,y) | x | y
#' whether taking lower triangle (default) or not
#' in this notation: x is index for row, y index for col
#'
#' @param hicMat a squared HiC matrix
#'
#' @return
#' @export
#'
#' @examples
convertForVtk <- function(hicMat) {
    hicMat_ <- as.matrix(hicMat)
    n <- dim(hicMat_)[2]
    ## hicMat_[lower.tri(hicMat_)] <- NA
    hicMat_ <- reshape2::melt(hicMat_, varnames = c('row', 'col'), na.rm = TRUE)
    colnames(hicMat_) <- c('x', 'y', 'f')
    ## get upTriId based on vtk format in python
    upTriId <- seq(0, n*(n-1), by = n)
    for (i in 1:(n-1)) {
        temp = upTriId + 1
        upTriId <- c(upTriId, tail(temp, n-i+1))
    }
    upTriId <- upTriId + 1
    return(hicMat_[upTriId,])
}

#' moving average
#'
#' @param x
#' @param order
#' @param center
#'
#' @return
#' @export
#'
#' @examples
ma <- function(x, order, center = TRUE)
{
    if (abs(order - round(order)) > 1e-08)
        stop("order must be an integer")
    if (order%%2 == 0 & center)
        w <- c(0.5, rep(1, order - 1), 0.5)/order
    else w <- rep(1, order)/order
    return(filter(x, w))
}

#' plot persistence curve
#'
#' @param pLevelGridList slot pLevelGrid from \code{treeHiCDataSet} object
#' @param path path from \code{treeHiCDataSet} object
#'
#' @return
#' @export
#'
#' @examples
plot_persistence_curve <- function(pLevelGridList, path) {
    pLevelGrid <- pLevelGridList[['pLevelGrid']]
    pLevel1 <- pLevelGrid[length(pLevelGrid)]
    nExtrema <- pLevelGridList[['nExtremaPLevel1']]
    dat_diff <- pLevelGridList[['dat_diff']]
    pLevel_diff <- pLevelGridList[['pLevel_diff']]
    dat <- pLevelGridList[['dat']]
    dfPlot <- list()
    ## ## persistence difference
    persCurvediff <- paste0('persistence vs. persistence difference: pLevel1 =', pLevel1, '; nExtremas = ', nExtrema)
    dfPlot[[1]] <- ggplot2::ggplot(dat_diff, aes(x = persistence, y = pers_diff)) +
        geom_point() +
        geom_hline(yintercept = pLevel_diff, col = 'blue') +
        geom_vline(xintercept = pLevel1, col = 'blue') + ggtitle(persCurvediff)
    dfPlot[[2]] <- ggplot2::ggplot(dat_diff, aes(x = persistence, y = pers_diff)) +
        geom_point() + coord_trans(x = "log10") +
        geom_hline(yintercept = pLevel_diff, col = 'blue') +
        geom_vline(xintercept = pLevel1, col = 'blue') + ggtitle("logx scale")
    ### persistence vs nExtremas
    persCurve <- paste0('log-log scale: pLevel1 =', pLevel1, '; nExtremas = ', nExtrema)
    dfPlot[[3]] <- ggplot2::qplot(y = nExtremas, x=persistence, data= dat, geom = "line", log = "xy", main = persCurve) +
        geom_vline(xintercept = pLevel1, col = 'blue')
    dfPlot[[4]] <- ggplot2::qplot(y = nExtremas, x=persistence, data= dat, geom = "point", log = "xy", main = persCurve) +
        geom_vline(xintercept = pLevel1, col = 'blue') +
        geom_hline(yintercept = nExtrema, col = 'blue')
    dfPlot[[5]] <- ggplot2::qplot(y = nExtremas, x=persistence, data= dat, geom = "line", log = "xy", main = persCurve) +
        geom_vline(xintercept = pLevel1, col = 'blue') +
        geom_vline(xintercept = pLevelGrid[1:(length(pLevelGrid)-1)], col = 'green')
    hist_title <- paste0('hist of log persistence with pLevel1 = ', pLevel1 )
    dfPlot[[6]] <-ggplot2::ggplot(dat, aes(log10(persistence))) + geom_density(adjust = 1/10) +
        geom_vline(xintercept = log10(pLevel1), col = 'blue') + ggtitle(hist_title)
    pdf(file = paste0(path,'temp/pLevelGrid.pdf'))
    for (i in 1:length(dfPlot)) {
        print(dfPlot[[i]])
    }
    hist(log10(dat$persistence), breaks = 300, xlim = c(-7,0), main= hist_title)
    abline(v= log10(pLevel1) )
    dev.off()
}

#' get cell ID (or vertexID vtk pipeline) for pLevel1 and pLevel2 with inputs
#' from Topological Toolkit (ttk)
#'
#' @param pathThresholdPLevelList list of csv files for ThresholdPLevel data
#' @param pathMSComplexPLevel1 MSComplex partitions from pLevel1
#' @param d_height data frame of 3 columns storing xy coords and height/difference values
#'
#' @return slot \code{testingTree} for object \code{treeHiCDataSet} object
#' @export
#'
#' @examples
create_testing_tree <- function(pathThresholdPLevelList, pathMSComplexPLevel1, d_height) {
    thresPLevelList <- list()
    fNames <- paste0('thresholdPLevel', 1:length(pathThresholdPLevelList))
    for (i in 1:length(pathThresholdPLevelList)) {
        thresPLevelList[[fNames[i]]] <- data.table::data.table(readr::read_csv(pathThresholdPLevelList[[i]],
                            col_types = cols(`Coordinates:0` = col_skip(), `Coordinates:1` = col_skip(),
                                             `Coordinates:2` = col_skip(), `Points:0` = col_skip(), `Points:1` = col_skip(),
                                             `Points:2` = col_skip(), RegionSize = col_skip(), RegionSpan = col_skip(),
                                             RegionType = col_skip())))
    }

    ## get MS partition for pLevel1
    MSComplexPLevel1 <- readr::read_csv(pathMSComplexPLevel1,
                                col_types = cols(`Points:2` = col_skip(), heights = col_skip()))
    names(MSComplexPLevel1)[c(5,6)] <- c('x_vtk', 'y_vtk' )
    MSComplexPLevel1 <- data.table::data.table(MSComplexPLevel1)
    NodeTypeSelect <- c(0, 4)
    vertexIdPLevelList <- list()
    ## case PLevel1: treats as a separate case
    IdPLevel1 <-  thresPLevelList[[1]][NodeType %in% NodeTypeSelect]$VertexIdentifier + 1
    vertexIdPLevel1 <- d_height[IdPLevel1, c('x', 'y')] - 1
    vertexIdPLevel1$vertexId <- IdPLevel1
    names(vertexIdPLevel1) <- c('x_vtk', 'y_vtk', 'vertexId' )
    ## get partition labels for each points in pLevel1
    ## TODO: only works on DescendingManifold (aka max points) for now
    vertexIdPLevel1 <- merge(MSComplexPLevel1, vertexIdPLevel1)[, c('x_vtk', 'y_vtk', 'DescendingManifold', 'vertexId')]
    names(vertexIdPLevel1)[3] <- 'ManifoldId'
    ## this list contains ALL POINTS, NOT ONLY EXTREMAS
    vertexPartitions <- lapply(unique(vertexIdPLevel1$ManifoldId), function(l) {MSComplexPLevel1[DescendingManifold == l, c('x_vtk', 'y_vtk')]} )
    names(vertexPartitions) <- unique(vertexIdPLevel1$ManifoldId)
    vertexIdPLevelList[['pLevel1']] <- list('rootNode' = vertexIdPLevel1)
    ## case: pLevel > pLevel1
    if (length(thresPLevelList) > 1) {
        for (thresIndex in 2:length(thresPLevelList)) {
        ## get extrema ids for testing at pLevel(s)
        ## NodeType: 0:local min, 4: local max
        ## add 1 since VTK is counting from 0
        IdPLevel <-  thresPLevelList[[thresIndex]][NodeType %in% NodeTypeSelect]$VertexIdentifier + 1
        vertexIdPLevel <- d_height[IdPLevel, c('x', 'y')] - 1
        vertexIdPLevel$vertexId <- IdPLevel
        names(vertexIdPLevel) <- c('x_vtk', 'y_vtk', 'vertexId' )
        ## ## get partition labels for each points in pLevel(s)
        ## ## TODO: only works on DescendingManifold (aka max points) for now
        vertexIdPLevel <- lapply(1:length(vertexPartitions), function(p) {data.table::data.table(merge(vertexIdPLevel, vertexPartitions[[p]]))} )
        names(vertexIdPLevel) <- names(vertexPartitions)
        vertexIdPLevelList[[paste0('pLevel',thresIndex)]] <- vertexIdPLevel
    }
    }
    return(vertexIdPLevelList)
}


#' add p-values information to
#' \code{testingTree}
#'
#' @param testingTree a slot \code{testingTree} of object \code{treeHiCDataSet}
#' @param mat_pvals a squared matrix of (raw) p-values to perform testing
#'
#' @return
#' @export
#'
#' @examples
evalPvals <- function(testingTree, mat_pvals) {
    output_result <- function(mat_pvals, i.range, j.range, toPrint = TRUE) {
        res <- matrix(nrow = length(i.range), ncol = 3)
        for (i in 1:length(i.range) ) {
            if (toPrint) {
                print(paste('location:', i.range[i], ',',j.range[i]))
                print(paste('pvals = ', mat_pvals[i.range[i], j.range[i]]))
            }
            res[i, ] <- c(i.range[i], j.range[i], mat_pvals[i.range[i], j.range[i]])
        }
        res <- data.frame(res)
        ## transform to vtk coordinates (off by 1)
        res[, 1] <- res[, 1] - 1
        res[, 2] <- res[, 2] - 1
        names(res) <- c('x_vtk', 'y_vtk', 'pvalues')
        return(data.table::data.table(res))
    }
    treeDepth <- length(testingTree)
    if (treeDepth == 0) {
        stop("tree is Empty (treeDepth = 0)")
    }
    for (treePLevel in 1:treeDepth) {
        totalNodes <- length(testingTree[[treePLevel]])
        if (totalNodes > 0 ) {
            for (ManifoldId in 1:totalNodes) {
                tree_level <- testingTree[[treePLevel]][[ManifoldId]]
                if (dim(tree_level)[1] > 0) {
                    ms_pvals <- output_result(mat_pvals = mat_pvals,
                                             i.range = tree_level$x_vtk + 1,
                                             j.range = tree_level$y_vtk + 1 ,
                                             toPrint = FALSE)
                    ms_pvals[, `:=`(p.adj, p.adjust(pvalues, method = 'BH'))]
                    testingTree[[treePLevel]][[ManifoldId]] <- merge(tree_level, ms_pvals)
                }
            }
        }
    }
    return(testingTree)
}

#' Title
#'
#' @param testingTree a slot \code{testingTree} of object \code{treeHiCDataSet}
#' @param alpha FDR control at \code{alpha} level
#' @param use_adjusted_pvals Logical, default to TRUE if adjusted pvalues are used
#'
#' @return
#' @export
#'
#' @examples
get_diff_hic <- function(testingTree, alpha, use_adjusted_pvals = TRUE) {
    print(paste0('adjusted is ', use_adjusted_pvals))
    hic_diff_res <- data.frame(matrix(nrow=0,ncol=6))
    names(hic_diff_res) <- c('x_vtk' , 'y_vtk', 'ManifoldId', 'vertexId', 'pvalues', 'p.adj')
    ## this keeps track of testing result in the tree (0: do not reject; 1: reject)
    checkTree <- list()
    treeDepth <- length(testingTree)
    if (treeDepth == 0) {
        stop("tree is Empty (treeDepth = 0)")
    }
    ## initialize checkTree
    ManifoldId <- as.character(unique(testingTree[[1]][[1]]$ManifoldId))
    for (treePLevel in 1:treeDepth) {
        checkTree[[treePLevel]] <- data.table("ManifoldId" = ManifoldId, "test_result" = rep(0, length(ManifoldId)))
    }
    names(checkTree) <- paste0('pLevel',1:treeDepth)
    for (treePLevel in 1:treeDepth) {
        if (treePLevel == 1) {
            tree_level <- testingTree[[treePLevel]][[1]]
            if (dim(tree_level)[1] > 0) {
                res <- tree_level[p.adj <= alpha]
                if (dim(res)[1] > 0) {
                    nonNull <- unique(res$ManifoldId)
                    checkTree[[treePLevel]][ManifoldId %in% nonNull]$test_result <- rep(1, length(nonNull))
                    hic_diff_res <- rbind(hic_diff_res, res)
                }
            } else {
                stop("number of nodes at pLevel1 are zero. Nothing to test ")
            }
        } else { # pLevel > pLevel1
            parents <- as.character(checkTree[[treePLevel - 1]][test_result == 1]$ManifoldId)
            if (length(parents) > 0) {
                res <- list()
                for (p in parents) {
                    tree_level <- testingTree[[treePLevel]][[p]]
                    if (dim(tree_level)[1] > 0) {
                        res[[p]] <- tree_level[p.adj <= alpha]
                    } else { #number of children at the current pLevel is empty: label 1? (not 0)
                        ##TODO: should this partition stops when its children at the current pLevel empty?
                        #checkTree[[treePLevel]][ManifoldId == p]$test_result <- 1
                    }
                }
                # nonNull: non-empty at thresholdPLevel, AND non-empty at [p.adj <= alpha]
                nonNull <- names(res)[sapply(res, dim)[1,] > 0]
                nonNull <- as.numeric(unique(nonNull))
                checkTree[[treePLevel]][ManifoldId %in% nonNull]$test_result <- rep(1, length(nonNull))
                for (p in as.character(nonNull)) {
                    dat <- data.frame(res[[p]], 'ManifoldId' = rep( as.numeric(p), nrow(res[[p]])))
                    dat <- dat[,c(1,2,6,3,4,5)]
                    hic_diff_res <- rbind(hic_diff_res, dat)
                }
            } else {
                cat("number of parent nodes are zeros. Tree stopped.")
            }
        }
    }
    if (dim(hic_diff_res)[1] > 0) {
        hic_diff_res[, 'x_vtk'] <- hic_diff_res[, 'x_vtk'] + 1
        hic_diff_res[, 'y_vtk'] <- hic_diff_res[, 'y_vtk'] + 1
        names(hic_diff_res)[c(1,2)] <- c('x', 'y')
        hic_diff_res <- as.matrix(hic_diff_res)
    }
    return(list('hic_diff_res' = hic_diff_res, 'checkTree' = checkTree))
}


#### ttk python pipeline ####

#' generate persistence curve
#'
#' @param path path to a working folder
#'
#' @return
#' @export
#'
#' @examples
get_persistence_curve <- function(path) {
    print("create vtk file")
    ## create-vtk.pt
    pyPath <- system.file("python", "create-vtk.py", package = "TreeHiC")
    ## python create-vtk.py path-to-temp-folder
    command <- paste0("python ", pyPath, ' ',path, 'temp/')
    system(command)
    ## ttkPF.py
    print("need paraview pvpython: generate persistence curve")
    pyPath <- system.file("python", "ttkPD.py", package = "TreeHiC")
    ## python create-vtk.py path-to-temp-folder
    command <- paste0("pvpython ", pyPath, ' ', path, 'temp/')
    system(command)
}

#' generate Morse-Smale partitions
#'
#' @param path path to a working folder
#' @param pLevelGrid a list whose elements contain different persistence levels of
#' \code{testingTree}
#'
#' @return
#' @export
#'
#' @examples
get_partitions <- function(path, pLevelGrid) {
    ## ttkMSComplex.py
    print("need paraview pvpython: create MS partition labels for pLevel1")
    pyPath <- system.file("python", "ttkMSComplex.py", package = "TreeHiC")
    ## python ttkMSComplex.py path-to-temp-folder pLevel1 1
    command <- paste0("pvpython ", pyPath, ' ', path, 'temp/ ', pLevelGrid[length(pLevelGrid)], ' 1')
    system(command)
    ## ttkPDThreshold.py
    print("create grid of pLevel for the testing tree")
    pyPath <- system.file("python", "ttkPDThreshold.py", package = "TreeHiC")
    ## python ttkPDThreshold.py path-to-temp-folder pLevel1 1 thresholdPLevel1
    command <- paste0("pvpython ", pyPath, ' ', path, 'temp/ ', pLevelGrid[length(pLevelGrid)],
                     ' 1 thresholdPLevel1')
    system(command)
    thres_names <- paste0('thresholdPLevel', 1:length(pLevelGrid))
    if (length(pLevelGrid) > 1) {
        thres_grid_st <-  sapply(length(pLevelGrid):2, function(i) {
            paste0("pvpython ", pyPath, ' ', path, 'temp/ ',
                   pLevelGrid[i-1], ' ',pLevelGrid[i], ' ', thres_names[length(pLevelGrid) - i + 2])
        })
        for (st in thres_grid_st) {
            print(st)
            system(st)
        }
    }
}


#' @import data.table
NULL
#' @import methods
NULL
#' @import readr
NULL
