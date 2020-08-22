suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(dplyr)
  library(Seurat)
  library(pbapply)
  library(future)
  library(future.apply)
  library(stats)
  library(SingleCellExperiment)
})

# ------------- Pre-process
load("/Users/xu/Downloads/ipf_qc_genes_union_biogrid_intact.Robj")  # Load Seurat object into r (Name : ipf)

# --------------------------------------------------------------------------------
# --------------------------- Review When You Run for Whole Dataset -------------
# Subset Data for Local Running
# ipf <- subset(ipf, subset = nFeature_RNA > 800 & nFeature_RNA < 1000  & nCount_RNA < 2500) # 2626 samples
# ipf <- subset(ipf, features = sample(c(1:dim(ipf[["RNA"]]@data)[1]), 1000))  # Subset first 1000 features for testing : 500 features * 2626 samples

# Save and Export Small data
# # save(ipf, file = "/Users/xu/OneDrive/SingleCell/SEURAT/Wilcox/ipd_test.Robj")

# ---------- Filter one cell type  (Here we use B cell for testing)
# table(ipf$cell.type.ident)
# ipf <- subset(ipf, subset = cell.type.ident == "B")

# --------------------------------------------------------------------------------

# Convert Seurat object to SingleCellAssay object
tem <- as.SingleCellExperiment(ipf)
tem <- SceToSingleCellAssay(tem)

# Check the Col Data for the transformed data
# colData(tem)
# dim(assay(tem))
# assay(tem)[1:10, 1:10]

# ---------------- GLOBAL SCALING ------
# Recalculating the cellular detection rate (ngeneson)
cdr2 <- colSums(assay(tem)>0)
SummarizedExperiment::colData(tem)$cngeneson <- scale(cdr2) # Global Scaling

# Factor and Relevel Control group and IPF
cond <- factor(x = SummarizedExperiment::colData(tem)$disease.ident)
cond <- relevel(x = cond, ref = "Control")
SummarizedExperiment::colData(tem)$condition <- cond


# Construct the Loop
cell.type <- names(table(ipf@meta.data$cell.type.ident))
MAST_ZScore <- list()

for (i in 1:length(cell.type)){
  # Subset for One single Cell Type
  test <- subset(tem, subset = cell.type.ident == cell.type[i])
  
  # Fit the ZLM model
  zlmCond <- MAST::zlm(~ condition + cngeneson, sca = test)
  
  # Get the summary datatable
  summaryCond <- summary(object = zlmCond)
  summaryDt <- as.data.frame(summaryCond$datatable)
  
  MAST_ZScore[[i]] <- summaryDt %>% 
    dplyr::filter(component == "logFC" & contrast == "conditionIPF")
}

names(MAST_ZScore) <- cell.type





# ----------- Test Version if we scaling cell-specificly
# cells.1 <- WhichCells(object = ipf, idents = "B_Control") # 174
# cells.2 <- WhichCells(object = ipf, idents = "B_IPF")     # 465
# 
# MASTzScore <- function(
#   data.use,
#   cells.1,
#   cells.2,
#   verbose = TRUE
# ){
#   if (!PackageCheck('MAST', error = FALSE)) {
#     stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
#   }
#   
#   # Create
#   group.info <- data.frame(row.names = c(cells.1, cells.2))
#   
#   
#   
#   
#   
# }
# 
# group.info <- data.frame(row.names = c(cells.1, cells.2))
# group.info[cells.1, "group"] <- "Group1"
# group.info[cells.2, "group"] <- "Group2"
# group.info[, "group"] <- factor(x = group.info[, "group"])

