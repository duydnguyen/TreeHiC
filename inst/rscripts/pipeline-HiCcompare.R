## codes based on "~/Research/HiCDA/Simulations/R/pipeline_withRep.R"
## date: 2018-02-23
## update: 2018-02-{26}
rm(list =ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))

library(ggplot2)
library(TreeHiC)
## for diffHiC
library(diffHic)
library(InteractionSet)
library(edgeR)

source("utility_functions.R")

#### @ params
alpha <- 0.05
maxDepth <- 1
use_adjusted_pvals <- TRUE
path_to_sim <- "~/ResearchLocal/HiCDA/Simulations/DeriveData/HiCcompare/fc1.5_n1000/simDbList.RData"
sim <- 8
## working path
path <- "~/ResearchLocal/HiCDA/Simulations/DeriveData/HiCcompare/fc1.5_n1000/" #must end with /

system(paste0('mkdir ', path, 'temp'))
sample.ids <- c('mat1', 'mat2')
txt_files <- rep(1:2)
condition <- factor(c(rep(c('cond1', 'cond2'), 1 )))
coldata <- data.frame('SampleID' = sample.ids, 'Condition' = condition, 'files' = txt_files)

#### manually adjust pLevelGrid
pLevelSelect <- list()
## load simulation data
load(path_to_sim)
simDb <- simDbList[[sim]]
d <- simDb[['d']]
DA_truth <- simDb[['DA_truth']]
i.range <- DA_truth[,1]
j.range <- DA_truth[,2]
mat_pvals <- simDb[['mat_pvals']]

matInputs <- list('mat1' = simDb[['C1']],
                  'mat2' = simDb[['C2']])

hicDb <- new("treeHiCDataSet")
hicDb <- HiCDataSetFromMatrix(hicDb, contactMatrixList = matInputs,
                              colData = coldata, path = path)

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

print(paste("Sparsity: n=", dim(mat_pvals)[1]))
print(paste("number of interactions: nChanges=", dim(DA_truth)[1]))

#### manually adjust pLevelGrid
# pLevelSelect[[1]] <- c(10^{-0.9}, 10^{-0.8}) #sim1
# pLevelSelect[[2]] <- c(10^{-0.99}, 10^{-0.95}) #sim2
# pLevelSelect[[3]] <- c(10^{-1.2}, 10^{-1}) #sim3
# pLevelSelect[[4]] <- c(10^{-1.2}, 10^{-1}) #sim4
# pLevelSelect[[5]] <- c(10^{-1.3}, 10^{-1}) #sim5
# pLevelSelect[[6]] <- c(10^{-1.3}, 10^{-1.2}) #sim6
# pLevelSelect[[7]] <- c(10^{-1.15}, 10^{-1}) #sim7
pLevelSelect[[8]] <- c(10^{-1.05}, 10^{-1}) #sim8
hicDb@pLevelGrid[["pLevelGrid"]] <- pLevelSelect[[sim]]
TreeHiC::plot_persistence_curve(hicDb@pLevelGrid, path = hicDb@path)

#### generate partitions
TreeHiC::get_partitions(path = hicDb@path, pLevelGrid = hicDb@pLevelGrid[["pLevelGrid"]])

#### DA ####
hicDb <- hic_diff(hicDb, mat_pvals = mat_pvals,use_adjusted_pvals = use_adjusted_pvals, alpha = alpha)
fdr_table_tree <- eval_fdr_table(pvals = hicDb@hic_diff_result, hic_dim = nrow(mat_pvals),
                               i.range.truth = i.range, j.range.truth = j.range)
print(fdr_table_tree)

#### compare with only permutation ####
fdr_table_perm <- eval_fdr_table_permutation(mat_pvals = mat_pvals, i.range = i.range,
                           j.range = j.range, alpha = alpha,
                           use_adjusted_pvals = use_adjusted_pvals)
print(fdr_table_perm)
# raw p-vals
fdr_table_perm_raw <- eval_fdr_table_permutation(mat_pvals = mat_pvals, i.range = i.range,
                                             j.range = j.range, alpha = alpha,
                                             use_adjusted_pvals = FALSE)
print(fdr_table_perm_raw)

#### compare with diffHiC ####
results.r <- eval_diffHiC_result(simDb)
fdr_table_diffHiC <- eval_fdr_table_diffHiC(results.r, dim_hic = dim(simDb[['C1']])[1], i.range, j.range, alpha = alpha)
print(fdr_table_diffHiC)

#### try MS method with p-values from diffHiC: much better results
hicDb <- hic_diff(hicDb, mat_pvals = res_diffHiC_,use_adjusted_pvals = use_adjusted_pvals, alpha = alpha)
fdr_table_tree_diffhic <- eval_fdr_table(pvals = hicDb@hic_diff_result, hic_dim = nrow(mat_pvals),
                               i.range.truth = i.range, j.range.truth = j.range)

##### Output

fdr_tables <- rbind(fdr_table_perm_raw,
                    fdr_table_tree,
                    fdr_table_diffHiC,
                    fdr_table_tree_diffhic)
methods <- c("fdr_table_perm_raw",
                 "fdr_table_tree",
                 "fdr_table_diffHiC",
                 "fdr_table_tree_diffhic")
fdr_tables <- data.frame(fdr_tables, 'methods' = methods)
fdr_tables

save(fdr_tables, file = paste0(path,'fdr_tables_sim',sim,'.RData'))
