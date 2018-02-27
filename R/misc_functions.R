#' Eval fdr table for simulation
#'
#' @param pvals
#' @param hic_dim
#' @param i.range.truth
#' @param j.range.truth
#'
#' @return
#' @export
#'
#' @examples
eval_fdr_table <- function(pvals, hic_dim, i.range.truth, j.range.truth) {
    if (dim(pvals)[1] == 0) {
        stop('matrix pvals is EMPTY')
    }
    res <- gold_res <- matrix(0, nrow = hic_dim, ncol = hic_dim)
    for (k in 1:nrow(pvals)) {
        res[pvals[k, 1], pvals[k, 2] ] <- res[pvals[k, 2], pvals[k, 1] ] <- 1
    }
    for (k in 1:length(i.range.truth)) {
        gold_res[i.range.truth[k], j.range.truth[k]] <- gold_res[j.range.truth[k], i.range.truth[k] ] <-1
    }
    TP <- FP <- FN <- TN <- 0
    for (i in 1:nrow(res)) {
        for (j in i:ncol(res)) {
            if (res[i,j] == 1 & gold_res[i,j] == 1 ) {
                TP <- TP + 1
            }
            if (res[i,j] == 1 & gold_res[i,j] == 0 ) {
                FP <- FP + 1
            }
            if (res[i,j] == 0 & gold_res[i,j] == 1 ) {
                FN <- FN + 1
            }
            if (res[i,j] == 0 & gold_res[i,j] == 0 ) {
                TN <- TN + 1
            }
        }
    }
    #TP; FP; FN; TN
    ## number of rejections
    R <- TP + FP
    R <- ifelse(R>0, R, 0)
    ## eFDR
    eFDR <- ifelse(R>0, FP/R, 0)
    ## Sensitivity or True Positive Rate(TPR)
    TPR <- TP/(TP + FN)
    ## Specificity or True Negative Rate
    TNR <- TN/(TN+FP)
    ## Precision
    PPV <- ifelse(R>0, TP/(TP+FP), 0)

    return(data.frame('eFDR' = eFDR, 'Sensitivity' = TPR, 'Specificity' = TNR, 'Precision' = PPV))
}

#' calculate permutation p-values for each cell
#'
#' @param mat1 a HiC contact matrix for the first condition
#' @param mat2 a HiC contact matrix for the first condition
#'
#' @return a squared matrix storing pvalues from permutation
#' @export
#'
#' @examples
get_perm_pvals <- function(mat1, mat2) {
    mat1 <- HiCcompare::full2sparse(mat1)
    mat2 <- HiCcompare::full2sparse(mat2)
    hic_table <- HiCcompare::create.hic.table(mat1, mat2, chr = 'simChr', scale = FALSE)
    ## jointly normalize data for a single chromosome
    hic_table <- HiCcompare::hic_loess(hic_table, Plot = FALSE, Plot.smooth = FALSE)
    ## testing by permutation
    hic_table_raw <- HiCcompare::hic_diff(hic_table, Plot = FALSE, Plot.smooth = FALSE, diff.thresh = NA)
    print('+++DONE computing hic_table_raw')
    ## obtain the matrix of pvalues
    mat_pvals <- HiCcompare::sparse2full(hic_table_raw, hic.table = TRUE, column.name = "p.value")
    return(mat_pvals)
}

#' Title
#'
#' @param mat_pvals
#' @param i.range
#' @param j.range
#' @param alpha
#' @param use_adjusted_pvals
#'
#' @return
#' @export
#'
#' @examples
## eval_fdr_table_permutation <- function(mat_pvals, i.range, j.range, alpha, use_adjusted_pvals) {
##     ## get the entire p.values (of course, only upper triangular)
##     mat_pvals[lower.tri(mat_pvals, diag = FALSE)] <- NA
##     mat_pvals.a <- matrix(p.adjust(mat_pvals, method = 'BH'), nrow = nrow(mat_pvals),
##                          ncol = ncol(mat_pvals))
##     if (!use_adjusted_pvals) {
##         nonNull_id <- which(mat_pvals <= alpha, arr.ind = T)
##     } else {
##         nonNull_id <- which(mat_pvals.a <= alpha, arr.ind = T)
##     }
##     print(paste0('number of called interactions = ', dim(nonNull_id)[1]))
##     HiCcompare_res <- gold_res <- matrix(0, nrow = nrow(mat_pvals), ncol = ncol(mat_pvals))
##     for (k in 1:nrow(nonNull_id)) {
##         HiCcompare_res[nonNull_id[k, 1], nonNull_id[k, 2] ] <- HiCcompare_res[nonNull_id[k, 2], nonNull_id[k, 1] ] <-1
##     }
##     for (k in 1:length(i.range)) {
##         gold_res[i.range[k], j.range[k] ] <- gold_res[j.range[k], i.range[k] ] <-1
##     }
##     TP <- FP <- FN <- TN <- 0
##     for (i in 1:nrow(HiCcompare_res)) {
##         for (j in 1:ncol(HiCcompare_res)) {
##             if (HiCcompare_res[i,j] == 1 & gold_res[i,j] == 1 ) {
##                 TP <- TP + 1
##             }
##             if (HiCcompare_res[i,j] == 1 & gold_res[i,j] == 0 ) {
##                 FP <- FP + 1
##             }
##             if (HiCcompare_res[i,j] == 0 & gold_res[i,j] == 1 ) {
##                 FN <- FN + 1
##             }
##             if (HiCcompare_res[i,j] == 0 & gold_res[i,j] == 0 ) {
##                 TN <- TN + 1
##             }
##         }
##     }
##     TP <- TP/2; FP <- FP/2 ; FN <- FN/2; TN <- TN/2
##     ## number of rejections
##     R <- TP + FP
##     ## eFDR
##     eFDR <- FP/R
##     ## Sensitivity or True Positive Rate(TPR)
##     TPR <- TP/(TP + FN)
##     ## Specificity or True Negative Rate
##     TNR <- TN/(TN+FP)
##     ## Precision
##     PPV <- TP/(TP+FP)
##     fdr_table <- data.frame('eFDR' = eFDR, 'Sensitivity' = TPR, 'Specificity' = TNR, 'Precision' = PPV)
##     return(fdr_table)
## }
eval_fdr_table_permutation <- function(mat_pvals, i.range, j.range, alpha, use_adjusted_pvals) {
    ## take upper triangular (including the diagonal), and remove NAs
    temp_mat_pvals <- temp_mat_pvals.a <- HiCcompare::full2sparse(mat_pvals)
    if (!use_adjusted_pvals) {
        ## nonNull_id <- which(mat_pvals <= alpha, arr.ind = T)
        nonNull_id <- as.matrix(temp_mat_pvals[IF <= alpha])
    } else {
        temp_mat_pvals.a$IF <- p.adjust(temp_mat_pvals$IF, method = 'BH')
        ## nonNull_id <- which(mat_pvals.a <= alpha, arr.ind = T)
        nonNull_id <- as.matrix(temp_mat_pvals.a[IF <= alpha])
        print(dim(nonNull_id))
    }
    print(paste0('number of called interactions = ', dim(nonNull_id)[1]))
    if (dim(nonNull_id)[1] == 0) {
        return(data.frame('eFDR' = 0, 'Sensitivity' = 0, 'Specificity' = 1, 'Precision' = 0))
    }
    HiCcompare_res <- gold_res <- matrix(0, nrow = nrow(mat_pvals), ncol = ncol(mat_pvals))
    for (k in 1:nrow(nonNull_id)) {
        HiCcompare_res[nonNull_id[k, 1], nonNull_id[k, 2] ] <- HiCcompare_res[nonNull_id[k, 2], nonNull_id[k, 1] ] <-1
    }
    for (k in 1:length(i.range)) {
        gold_res[i.range[k], j.range[k] ] <- gold_res[j.range[k], i.range[k] ] <-1
    }
    TP <- FP <- FN <- TN <- 0
    for (i in 1:nrow(HiCcompare_res)) {
        for (j in 1:ncol(HiCcompare_res)) {
            if (HiCcompare_res[i,j] == 1 & gold_res[i,j] == 1 ) {
                TP <- TP + 1
            }
            if (HiCcompare_res[i,j] == 1 & gold_res[i,j] == 0 ) {
                FP <- FP + 1
            }
            if (HiCcompare_res[i,j] == 0 & gold_res[i,j] == 1 ) {
                FN <- FN + 1
            }
            if (HiCcompare_res[i,j] == 0 & gold_res[i,j] == 0 ) {
                TN <- TN + 1
            }
        }
    }
    TP <- TP/2; FP <- FP/2 ; FN <- FN/2; TN <- TN/2
    ## number of rejections
    R <- TP + FP
    R <- ifelse(R>0, R, 0)
    ## eFDR
    eFDR <- ifelse(R>0, FP/R, 0)
    ## Sensitivity or True Positive Rate(TPR)
    TPR <- TP/(TP + FN)
    ## Specificity or True Negative Rate
    TNR <- TN/(TN+FP)
    ## Precision
    PPV <- ifelse(R>0, TP/(TP+FP), 0)
    fdr_table <- data.frame('eFDR' = eFDR, 'Sensitivity' = TPR, 'Specificity' = TNR, 'Precision' = PPV)
    return(fdr_table)
}
