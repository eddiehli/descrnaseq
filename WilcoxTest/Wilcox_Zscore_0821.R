.libPaths("/gpfs/ysm/pi/zhao/hl732/R_pkgs")
library(dplyr)
library(Seurat)
library(pbapply)
library(future)
library(future.apply)
library(stats)
library(rlist)
library(MAST)
library(coin)

# Load data in R
load("/Users/xu/Downloads/ipf_qc_genes_union_biogrid_intact.Robj")
cell.type <- table(ipf@meta.data$cell.type.ident)
celltype.disease <- names(table(ipf$celltype.disease))

# Load test dataset
# load("/Users/xu/OneDrive/SingleCell/SEURAT/Wilcox/ipd_test.Robj")

# Wilcox_test ONLY for first gene "DMP1" of B cell
wilcoxonZ(x = ipf[["RNA"]]@data[1, WhichCells(object = ipf,idents = "B_Control")], 
          y = ipf[["RNA"]]@data[1, WhichCells(object = ipf,idents = "B_IPF")])


#### T-statistics
DiffTTest_sta <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE
) {
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  t_sta <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        wilcoxonZ(x = data.use[["RNA"]]@data[x, cells.1], 
                  y = data.use[["RNA"]]@data[x, cells.2])
      }
    )
  )
  to.return <- data.frame(t_sta,row.names = rownames(x = data.use))
  return(to.return)
}

results_tstat <- list()
for (i in 1:length(cell.type)){
  test_tstat <- DiffTTest_sta(ipf,
                              cells.1 = WhichCells(object = ipf, idents = celltype.disease[2*i-1]),
                              cells.2 = WhichCells(object = ipf, idents = celltype.disease[2*i]))
  results_tstat[[i]] <- test_tstat
}
names(results_tstat) = names(cell.type)
tsta_matrix <- matrix(unlist(results_tstat), ncol = 18 ,byrow = F)
colnames(tsta_matrix) <- list.names(results_tstat)
row.names(tsta_matrix) <- row.names(results_tstat[[1]])
results_tstat[["Matrix"]] <- tsta_matrix

# -------- Calculate p_value according to Z score
#  
#  2 * pnorm(-abs(z))



# Export and Save
saveRDS(results_tstat, file = "./data/results_tstat.rds")





# ----------------- Test -------------
# results_tstat <- list()
# for (i in 1:2){
#   test_tstat <- DiffTTest_sta(ipf,
#                               cells.1 = WhichCells(object = ipf, idents = celltype.disease[2*i-1]),
#                               cells.2 = WhichCells(object = ipf, idents = celltype.disease[2*i]))
#   results_tstat[[i]] <- test_tstat
# }



