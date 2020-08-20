.libPaths("/gpfs/ysm/pi/zhao/hl732/R_pkgs")
library(dplyr)
library(Seurat)
library(pbapply)
library(future)
library(future.apply)
library(stats)
library(rlist)
library(MAST)

load("/Users/xu/Downloads/ipf_qc_genes_union_biogrid_intact.Robj")
cell.type <- table(ipf@meta.data$cell.type.ident)
celltype.disease <- names(table(ipf$celltype.disease))


# T-TEST ONLY for first gene "DMP1" of B cell
wilcox.test(x = ipf[["RNA"]]@data[1, WhichCells(object = ipf,idents = "B_Control")], 
            y = ipf[["RNA"]]@data[1, WhichCells(object = ipf,idents = "B_IPF")])$p.value


#### ------------- p-values
DiffTTest_pval <- function(
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
  p_val <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        wilcox.test(x = data.use[["RNA"]]@data[x, cells.1], y = data.use[["RNA"]]@data[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(p_val,row.names = rownames(x = data.use))
  return(to.return)
}
results_pval <- list()
for (i in 1:length(cell.type)){
  test_pval <- DiffTTest_pval(ipf,
                              cells.1 = WhichCells(object = ipf, idents = celltype.disease[2*i-1]),
                              cells.2 = WhichCells(object = ipf, idents = celltype.disease[2*i]))
  results_pval[[i]] <- test_pval
}
names(results_pval) = names(cell.type)
pval_matrix <- matrix(unlist(results_pval), ncol = 18 ,byrow = F)
colnames(pval_matrix) <- list.names(results_pval)
row.names(pval_matrix) <- row.names(results_pval[[1]])
results_pval[["Matrix"]] <- pval_matrix

# ---Export and Save
saveRDS(results_pval, file = "./data/results_pval.rds")




