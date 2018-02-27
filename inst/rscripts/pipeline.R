rm(list =ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))

library(ggplot2)
library(TreeHiC)
# library(devtools)
# document()
# load_all()

## params
alpha <- 0.05
maxDepth <- 1
use_adjusted_pvals <- TRUE
# path_to_sim <- "~/ResearchLocal/HiCDA/Simulations/DeriveData/2018-02-16_simulation/n1000_100_fc2/simDb3.RData"
path_to_sim <- "~/ResearchLocal/HiCDA/Simulations/DeriveData/2018-02-17_simulation/simDb2.RData"

## working path
path <- "~/ResearchLocal/HiCDA/test-TreeHiC/" #must end with /
system(paste0('mkdir ', path, 'temp'))

sample.ids <- c('rep1', 'rep2')
txt_files <- rep(1:2)
condition <- factor(c(rep(c('cond1', 'cond2'), 1 )))
coldata <- data.frame('SampleID' = sample.ids, 'Condition' = condition, 'files' = txt_files)

## load simulation data
load(path_to_sim)
d <- simDb[['d']]
DA_truth <- simDb[['DA']]
i.range <- DA_truth[,1]
j.range <- DA_truth[,2]

matInputs <- list('rep1' = simDb[['rep1']],
                  'rep2' = simDb[['rep2']])

hicDb <- new("treeHiCDataSet")
hicDb <- HiCDataSetFromMatrix(hicDb, contactMatrixList = matInputs,
                              colData = coldata, path = path)

# create mat_pvals from HiCCompare (not belong to TreeHiC package!)
mat_pvals <- TreeHiC::get_perm_pvals(simDb[['rep1']], simDb[['rep2']])

#### eval height function
hicDb <- evalDiffMat(hicDb)
## write.csv
write.csv(hicDb@d_height[,'f'], file = paste0(hicDb@path,"temp/heights-scalar.csv"), row.names = FALSE, quote=FALSE)

#### generate persistence curve
TreeHiC::get_persistence_curve(path = hicDb@path)
#### selectPLevelGrid
hicDb <- selectPLevelGrid(hicDb, maxDepth = maxDepth)
## plot
TreeHiC::plot_persistence_curve(hicDb@pLevelGrid, path = hicDb@path)

# ## manually adjust pLevelGrid
# hicDb@pLevelGrid[["pLevelGrid"]] <- c(10^{-3}, 10^{-2})

#### generate partitions
TreeHiC::get_partitions(path = hicDb@path, pLevelGrid = hicDb@pLevelGrid[["pLevelGrid"]])

#### DA ####
hicDb <- hic_diff(hicDb, mat_pvals = mat_pvals,use_adjusted_pvals = use_adjusted_pvals, alpha = alpha)
fdr_table_ms <- eval_fdr_table(pvals = hicDb@hic_diff_result, hic_dim = nrow(mat_pvals),
                               i.range.truth = i.range, j.range.truth = j.range)
print(fdr_table_ms)


#### compare with only permutation ####
fdr_table_perm <- eval_fdr_table_permutation(mat_pvals = mat_pvals, i.range = i.range,
                           j.range = j.range, alpha = alpha, use_adjusted_pvals = use_adjusted_pvals)
print(fdr_table_perm)

#### compare with diffHiC ####
# https://support.bioconductor.org/p/90184/
library(diffHic)
library(InteractionSet)
library(edgeR)

gr <- GRanges(seqnames = "chrSim",
              ranges = IRanges(start = 1:dim(simDb[['rep1']])[1], width = 2))
cm1 <- InteractionSet::ContactMatrix(simDb[['rep1']], gr, gr)
cm2 <- InteractionSet::ContactMatrix(simDb[['rep2']], gr, gr)
## filtering 0?
# to.keep <- as.matrix(cm1)!=0 | as.matrix(cm2)!=0 # use this to exclude 0
to.keep <- as.matrix(cm1)!=0 |  as.matrix(cm2)!=0 | as.matrix(cm1)==0 | as.matrix(cm2)==0
iset1 <- deflate(cm1, extract=to.keep)
iset2 <- deflate(cm2, extract=to.keep)

data <- cbind(iset1, iset1,iset2, iset2)
interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")

## construct the design matrix
design <- model.matrix(~ factor(c(1,1,2,2)))

## filtering out uninteresting bin pairs
# keep <- aveLogCPM(asDGEList(data)) > 0
# data <- data[keep,]

## normalizing counts between libraries
data <- normalize(data, type="loess")
y <- asDGEList(data)
## modelling biological variability
y <- estimateDisp(y, design = design)
fit <- glmQLFit(y, design = design)
## 6. testing for significant differences between groups
result <- glmQLFTest(fit)
#clusters <- diClusters(data, result$table, target=0.05, cluster.args=list(tol=1))
#topTags(result)
#rowData(data) <- cbind(rowData(data), result$table)
adj.p <- p.adjust(result$table$PValue, method="BH")
sum(adj.p <= alpha)

useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
inter.frame <- as.data.frame(interactions(data))[,useful.cols]
results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
o.r <- order(results.r$PValue)
#write.table(results.r[o.r,], file="~/Downloads/2018-02-14/binpairs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
res_diffHiC <- data.table::data.table(results.r[o.r,])

names(res_diffHiC)[10] <- "p.value" # raw p.value
res_diffHiC_ <- HiCcompare::sparse2full(res_diffHiC, hic.table = TRUE, column.name = "p.value")


#### eval FDR of diffHiC
idDA <- which(adj.p <= alpha)
diffHiC_res <- gold_res <- matrix(0, nrow = nrow(mat_pvals), ncol = ncol(mat_pvals))
for (k in 1:length(idDA)) {
    diffHiC_res[results.r[idDA[k],'start1'], results.r[idDA[k],'start2'] ] <- diffHiC_res[results.r[idDA[k],'start2'], results.r[idDA[k],'start1'] ] <-1
}
for (k in 1:length(i.range)) {
    gold_res[i.range[k], j.range[k] ] <- gold_res[j.range[k], i.range[k] ] <-1
}
TP <- FP <- FN <- TN <- 0
for (i in 1:nrow(diffHiC_res)) {
    for (j in 1:ncol(diffHiC_res)) {
        if (diffHiC_res[i,j] == 1 & gold_res[i,j] == 1 ) {
            TP <- TP + 1
        }
        if (diffHiC_res[i,j] == 1 & gold_res[i,j] == 0 ) {
            FP <- FP + 1
        }
        if (diffHiC_res[i,j] == 0 & gold_res[i,j] == 1 ) {
            FN <- FN + 1
        }
        if (diffHiC_res[i,j] == 0 & gold_res[i,j] == 0 ) {
            TN <- TN + 1
        }
    }
}
TP <- TP/2; FP <- FP/2 ; FN <- FN/2; TN <- TN/2
## number of rejections
R <- TP + FP
## eFDR
eFDR <- FP/R
## Sensitivity or True Positive Rate(TPR)
TPR <- TP/(TP + FN)
## Specificity or True Negative Rate
TNR <- TN/(TN+FP)
## Precision
PPV <- TP/(TP+FP)
fdr_table_diffHiC <- data.frame('eFDR' = eFDR, 'Sensitivity' = TPR, 'Specificity' = TNR, 'Precision' = PPV)
print(fdr_table_diffHiC)

#### try MS method with p-values from diffHiC: much better results
hicDb <- hic_diff(hicDb, mat_pvals = res_diffHiC_,use_adjusted_pvals = use_adjusted_pvals, alpha = alpha)
fdr_table_ms_diffhic <- eval_fdr_table(pvals = hicDb@hic_diff_result, hic_dim = nrow(mat_pvals),
                               i.range.truth = i.range, j.range.truth = j.range)

##### Output
print('order: perm, ms_perm, diffHiC, ms_diffHiC')
fdr_table_perm
fdr_table_ms
fdr_table_diffHiC
fdr_table_ms_diffhic
