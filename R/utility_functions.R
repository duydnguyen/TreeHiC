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
    hicMat_ <- reshape2::melt(hicMat_, varnames = c('row', 'col'), na.rm = TRUE)
    colnames(hicMat_) <- c('x', 'y', 'f')
    point_id_upper_tri <- upper.tri(hicMat, diag = TRUE)
    offDiagId <- seq(2, n^2-n, by = n+1)
    point_id_upper_tri[offDiagId] <- rep(TRUE, n-1)
    upTriId <- which(point_id_upper_tri == TRUE)
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
    return(stats::filter(x, w))
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
    ## for (i in 1:length(dfPlot)) {
    for (i in 6) {
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
#' @param batch_mode logical; TRUE if running on large HiC matrix. Default to FALSE
#' @param path_vtk_coords path to vtk-coords.csv; use this only if \code{batch_mode} is TRUE.
#'
#' @return slot \code{testingTree} for object \code{treeHiCDataSet} object
#' @export
#'
#' @examples
create_testing_tree <- function(pathThresholdPLevelList, pathMSComplexPLevel1, d_height,
                               batch_mode = FALSE, path_vtk_coords = "") {
    thresPLevelList <- list()
    ## extract (x_vtk, y_vtk) from vtk-coords.csv
    if (batch_mode) {
        vtk_coords <- readr::read_csv(path_vtk_coords, col_names = FALSE)
    }
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
    if (!batch_mode) {
        vertexIdPLevel1 <- d_height[IdPLevel1, c('x', 'y')] - 1
        vertexIdPLevel1$vertexId <- IdPLevel1
    } else {
        ## vtk_coords <- readr::read_csv(path_vtk_coords, col_names = FALSE)
        vertexIdPLevel1 <- vtk_coords[IdPLevel1, ]
        vertexIdPLevel1$vertexId <- IdPLevel1
    }
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
            if (!batch_mode) {
                vertexIdPLevel <- d_height[IdPLevel, c('x', 'y')] - 1
                vertexIdPLevel$vertexId <- IdPLevel
            } else {
                ## vtk_coords <- readr::read_csv(path_vtk_coords, col_names = FALSE)
                vertexIdPLevel <- vtk_coords[IdPLevel, ]
                vertexIdPLevel$vertexId <- IdPLevel
            }
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
#' @param batch_mode logical; TRUE if running on a large HiC matrix. Default to FALSE
#'
#' @return
#' @export
#'
#' @examples
get_diff_hic <- function(testingTree, alpha,
                        use_adjusted_pvals = TRUE, batch_mode = FALSE) {
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
                if (use_adjusted_pvals) {
                    res <- tree_level[p.adj <= alpha]
                } else {
                    res <- tree_level[pvalues <= alpha]
                }
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
                        if (use_adjusted_pvals) {
                            res[[p]] <- tree_level[p.adj <= alpha]
                        } else {
                            res[[p]] <- tree_level[pvalues <= alpha]
                        }
                    } else { #number of children at the current pLevel is empty: label 1? (not 0)
                        #### *TODO: should this partition stops when its children at the current pLevel empty?
                        ##checkTree[[treePLevel]][ManifoldId == p]$test_result <- 1
                        res[[p]] <- data.frame()
                    }
                }
                ## nonNull: non-empty at thresholdPLevel, AND non-empty at [p.adj <= alpha]
                nonNull <- names(res)[sapply(res, dim)[1,] > 0]
                if (!batch_mode) {
                    nonNull <- as.numeric(unique(nonNull))
                }
                checkTree[[treePLevel]][ManifoldId %in% nonNull]$test_result <- rep(1, length(nonNull))
                for (p in as.character(nonNull)) {
                    if (!batch_mode) {
                        dat <- data.frame(res[[p]], 'ManifoldId' = rep( as.numeric(p), nrow(res[[p]])))
                    } else {
                        dat <- data.frame(res[[p]], 'ManifoldId' = rep(p, nrow(res[[p]])))
                    }
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
#' @param eval_MSComplex logical to indicate if MS partition is evaluated.
#'
#' @return
#' @export
#'
#' @examples
get_partitions <- function(path, pLevelGrid, eval_MSComplex = TRUE) {
    if (eval_MSComplex) {
        ## ttkMSComplex.py
        print("need paraview pvpython: create MS partition labels for pLevel1")
        pyPath <- system.file("python", "ttkMSComplex.py", package = "TreeHiC")
        ## python ttkMSComplex.py path-to-temp-folder pLevel1 1
        command <- paste0("pvpython ", pyPath, ' ', path, 'temp/ ', pLevelGrid[length(pLevelGrid)], ' 1')
        system(command)
    }
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

#'  set up data for running batch version
#'
#' @param hicDb a \code{treeHiCDataSet} object
#' @param batchSize number of cells in batches
#' @param nBatches number of batches
#'
#' @return batches of csv files and hicDb objects
#' @export
#' @examples
create_batches <- function(hicDb, batchSize = 10^6, nBatches = 5) {
    hicDbList <- list()
    points_id_select_all <- rownames(hicDb@d_height)
    max_grid <- dim(hicDb@d_height)[1] / batchSize
    max_grid <- ifelse(max_grid > 1, floor(max_grid), 0)
    if (max_grid > 0) {
        grid_range <- c(0, seq(1, max_grid, by= 1) * batchSize, dim(hicDb@d_height)[1] )
    } else {
        stop("Dividing grids is not required since batchSize > number of cells tested.")
    }
    df_region_select <- data.frame('i' = c(), "min"=c(), "max"=c())
    for (i in 2:length(grid_range)) {
        batch_range <-  grid_range[i-1]:grid_range[i]
        points_id_select <- points_id_select_all[batch_range]
        df_region_select <- rbind(df_region_select,
                                 data.frame('i' = i,
                                            'min'= min(hicDb@d_height[batch_range,'f']),
                                            'max'= max(hicDb@d_height[batch_range,'f'])))
    }
    ## exclude regions whose cells are zeros
    df_region_select <- data.table::data.table(df_region_select)
    df_region_select <- df_region_select[min != max]
    grid_index_select <- df_region_select$i
#### write to csv for each batches
    if (nBatches == 1) {
        hicDb_temp <- new("treeHiCDataSet")
        df_csv <- data.frame('f' = c(), "points_id_select" = c())
        for (i in grid_index_select) {
            print(i)
            start <- grid_range[i-1] + 1
            end <- grid_range[i]
            df_csv <- rbind(df_csv,
                           data.frame('f' = hicDb@d_height[start:end, 'f'],
                                      'points_id_select' = points_id_select_all[start:end]))
           }
        hicDbList[[batch]] <- hicDb_temp
        write.csv(df_csv,
                  file = paste0(hicDb@path,"temp/heights-scalar",".csv"),
                  row.names = FALSE, quote=FALSE)
    } else { ## run multiple batches (nBatches > 1)
        print("run multiple batches")
        sizeBatch <- floor(length(grid_index_select) / nBatches)
        batchs_index <- list()
        id_seq <- ceiling(seq(1, length(grid_index_select), length.out = nBatches+1))
        for (batch in 1:nBatches) {
            start <- id_seq[batch]
            end <- id_seq[batch+1] - 1
            if (batch == nBatches) {
                end <- id_seq[batch+1]
            }
            batchs_index[[batch]] <- grid_index_select[start:end]
        }
        print(batchs_index)
        for (batch in 1:nBatches) {
            print(paste0('batch = ', batch))
            system(paste0('mkdir ', hicDb@path, 'temp/batch', batch))
            system(paste0('mkdir ', hicDb@path, 'temp/batch', batch,'/temp'))
            hicDb_temp <- new("treeHiCDataSet")
            batch_range_grid <- batchs_index[[batch]]
            df_csv <- data.frame('f' = c(), "points_id_select" = c())
            for (i in batch_range_grid) {
                start <- grid_range[i-1] + 1
                end <- grid_range[i]
                df_csv <- rbind(df_csv,
                               data.frame('f' = hicDb@d_height[start:end, 'f'],
                                          'points_id_select' = points_id_select_all[start:end]))
               }
            hicDbList[[batch]] <- hicDb_temp
            write.csv(df_csv,
                      file = paste0(hicDb@path,"temp/batch",batch,"/temp/","heights-scalar",".csv"),
                      row.names = FALSE, quote=FALSE)
        }
    }
    return(hicDbList)
}

#' combine hic_diff from all batches
#'
#' @param hicDbList a list of \code{treeHiCDataSet} objects resulting from batch mode
#'
#' @return a list (\code{testingTree}) which is a combined result from individual batch's \code{testingTree}.
#' @export
#'
#' @examples
 combine_testingTree <- function(hicDbList) {
    testingTree <- list()
    nBatches <- length(hicDbList)
    treeDepth <- length(hicDbList[[1]]@testingTree)
    if (treeDepth == 0) {
        stop("tree is Empty (treeDepth = 0)")
    }
    for (treePLevel in 1:treeDepth) {
        if (treePLevel == 1) {
            tree_level <- data.frame('x_vtk'=c(), 'y_vtk'=c(), 'ManifoldId'=c(),
                                    'vertexId'=c(), 'pvalues'=c(), 'p.adj'=c())
            for (batch in 1:nBatches) {
                tree_level_temp <- hicDbList[[batch]]@testingTree[[treePLevel]][[1]]
                ## relabel ManifoldId
                temp_lab <- paste0(rep(batch, length(tree_level_temp$ManifoldId)),'_',
                                  tree_level_temp$ManifoldId)
                tree_level_temp$ManifoldId <- temp_lab
                tree_level <- rbind(tree_level, tree_level_temp)
            }
            tree_level$p.adj <- NULL
            tree_level[, `:=`(p.adj, p.adjust(pvalues, method = 'BH'))]
            testingTree[[treePLevel]] <- list('rootNode' = tree_level)
        } else { ## treePLevel > 1
            tree_level <- list()
            for (batch in 1:nBatches) {
                tree_level_temp <- hicDbList[[batch]]@testingTree[[treePLevel]]
                ## relabel ManifoldId
                temp_lab <- paste0(rep(batch, length(names(tree_level_temp))),'_',
                                  names(tree_level_temp))
                names(tree_level_temp) <- temp_lab
                tree_level <- c(tree_level, tree_level_temp)
            }
            testingTree[[treePLevel]] <- tree_level
        }
    }
    names(testingTree) <- paste0('pLevel', 1:length(testingTree))
    return(testingTree)
}

#' get cell ID (or vertexID vtk pipeline) for pLevel1 and pLevel2 with inputs
#' from Topological Toolkit (ttk). This is the simplified version of hic_diff()
#' where the input of MSComplex is not required (i.e.,\code{eval_MSComplex = F}).
#'
#' @param object : a \code{treeHiCDataSet} object
#' @param mat_pvals : a squared matrix of (raw) p-values to perform testing
#' @param use_adjusted_pvals : Logical, should adjusted p-values be used?
#' @param alpha : FDR control at \code{alpha} level
#' @param batch_mode logical; TRUE if running on large HiC matrix. Default to FALSE
#'
#' @return add additional slots \code{testingTree, hic_diff_result}
#' @export
#'
#' @examples
hic_diff_simplified <- function(hicDb, mat_pvals, use_adjusted_pvals = TRUE,
                               alpha, batch_mode = FALSE) {
    ## create_testing_tree_simplified()
    create_testing_tree_simplified <- function(pathThresholdPLevelList,
                                              batch_mode = FALSE, path_vtk_coords = "") {
        thresPLevelList <- list()
        ## extract (x_vtk, y_vtk) from vtk-coords.csv
        if (batch_mode) {
            vtk_coords <- readr::read_csv(path_vtk_coords, col_names = FALSE)
        }
        fNames <- paste0('thresholdPLevel', 1:length(pathThresholdPLevelList))
        for (i in 1:length(pathThresholdPLevelList)) {
            thresPLevelList[[fNames[i]]] <- data.table::data.table(readr::read_csv(pathThresholdPLevelList[[i]],
                            col_types = cols(`Coordinates:0` = col_skip(), `Coordinates:1` = col_skip(),
                                             `Coordinates:2` = col_skip(), `Points:0` = col_skip(), `Points:1` = col_skip(),
                                             `Points:2` = col_skip(), RegionSize = col_skip(), RegionSpan = col_skip(),
                                             RegionType = col_skip())))
        }
        NodeTypeSelect <- c(0, 4)
        vertexIdPLevelList <- list()
        ## case PLevel1: treats as a separate case
        IdPLevel1 <-  thresPLevelList[[1]][NodeType %in% NodeTypeSelect]$VertexIdentifier + 1
        if (!batch_mode) {
            vertexIdPLevel1 <- d_height[IdPLevel1, c('x', 'y')] - 1
            vertexIdPLevel1$vertexId <- IdPLevel1
        } else {
            vertexIdPLevel1 <- vtk_coords[IdPLevel1, ]
            vertexIdPLevel1$vertexId <- IdPLevel1
        }
        names(vertexIdPLevel1) <- c('x_vtk', 'y_vtk', 'vertexId' )
        vertexIdPLevelList[['pLevel1']] <- list('rootNode' = vertexIdPLevel1)
        return(vertexIdPLevelList)
    }
    ## get_diff_hic_simplified()
    get_diff_hic_simplified <- function(testingTree, alpha,
                                       use_adjusted_pvals = TRUE, batch_mode = FALSE) {
        print(paste0('adjusted is ', use_adjusted_pvals))
        hic_diff_res <- data.frame(matrix(nrow=0,ncol=5))
        ## names(hic_diff_res) <- c('x_vtk' , 'y_vtk', 'ManifoldId', 'vertexId', 'pvalues', 'p.adj')
        names(hic_diff_res) <- c('x_vtk' , 'y_vtk', 'vertexId', 'pvalues', 'p.adj')
        ## this keeps track of testing result in the tree (0: do not reject; 1: reject)
        checkTree <- list()
        treeDepth <- length(testingTree)
        if (treeDepth == 0) {
            stop("tree is Empty (treeDepth = 0)")
        }
        for (treePLevel in 1:treeDepth) {
            if (treePLevel == 1) {
                tree_level <- data.table::data.table(testingTree[[treePLevel]][[1]])
                if (dim(tree_level)[1] > 0) {
                    if (use_adjusted_pvals) {
                        res <- tree_level[p.adj <= alpha]
                    } else {
                        res <- tree_level[pvalues <= alpha]
                    }
                    if (dim(res)[1] > 0) {
                        hic_diff_res <- rbind(hic_diff_res, res)
                    }
                } else {
                    stop("number of nodes at pLevel1 are zero. Nothing to test ")
                }
            } else { # pLevel > pLevel1
            }
        }
        if (dim(hic_diff_res)[1] > 0) {
            hic_diff_res[, 'x_vtk'] <- hic_diff_res[, 'x_vtk'] + 1
            hic_diff_res[, 'y_vtk'] <- hic_diff_res[, 'y_vtk'] + 1
            names(hic_diff_res)[c(1,2)] <- c('x', 'y')
            hic_diff_res <- as.matrix(hic_diff_res)
        }
        return(list('hic_diff_res' = hic_diff_res))
    }
    ## MAIN ##
    path_temp <- paste0(hicDb@path,"temp/")
    pLevelGrid <- hicDb@pLevelGrid[['pLevelGrid']]
    thres_names <- paste0('thresholdPLevel', 1:length(pLevelGrid))
    path_vtk_coords <- paste0(path_temp, "vtk-coords.csv")
    pathThresholdPLevelList <- list()
    for (i in 1:length(pLevelGrid)) {
        pathThresholdPLevelList[[thres_names[i]]] <- paste0(path_temp, paste0(thres_names[i],'.csv'))
    }
    ## create testing tree
    testingTree <- create_testing_tree_simplified(pathThresholdPLevelList = pathThresholdPLevelList,
                                               batch_mode = batch_mode, path_vtk_coords = path_vtk_coords)
    testingTree <- TreeHiC::evalPvals(testingTree = testingTree, mat_pvals = mat_pvals)
    ## perform differential testing
    diffResultList <- get_diff_hic_simplified(testingTree = testingTree, alpha = alpha, use_adjusted_pvals = use_adjusted_pvals)
    if (dim(diffResultList[['hic_diff_res']])[1] > 0) {
        hicDb@hic_diff_result <- as.matrix(diffResultList[['hic_diff_res']])
    }
    ## return hicDb
    hicDb@testingTree <- testingTree
    hicDb
}


#' @import data.table
NULL
#' @import methods
NULL
#' @import readr
NULL
