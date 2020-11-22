library(dplyr)
library(GenomicRanges)
library(plyranges)
library(WGCNA)
library(readr)
library(readxl)
library(fastcluster)
library(flashClust)
library(ComplexHeatmap)
library(data.table)
library(tidyverse)
library(circlize)

# ----------------------- Read Data ------------
# Load in the Result Data
load("/Users/xu/OneDrive/SingleCell/Heatmap for Region/data.ratio.normalized.snpFiltered.invariantFiltered.gr.RData")
load("/Users/xu/OneDrive/SingleCell/Heatmap for Region/data.annotation.gr.RData")
group <- read_excel("/Users/xu/OneDrive/SingleCell/CPG/MetadataFinal041720.xlsx")
# Address
# which(data.ratio.gr@elementMetadata %>% colnames %in% meta_dat$SampleID == 0)
colnames(data.ratio.gr@elementMetadata)[77] = "L5236_hipA"
colnames(data.ratio.gr@elementMetadata)[188] = "L5626_amygA"
# sum(!data.ratio.gr@elementMetadata %>% colnames %in% meta_dat$SampleID)

# Region metadata$RegionSc MNA
# Gene: ELFN1, LHX2, MGMT, PTDSS2, ADD3, CUEDC2, HARS2, PI4KB

# 102 Sample but only 99 in raw data
mata_sub <- group %>%
  filter(RegionSci == "MNA") %>%
  filter(DxGroup1 %in% c("PTSD", "Control"))


gene_8 <- data.ratio.gr %>% 
  select(any_of(mata_sub$SampleID))

# Subset annotation using GENE names
gene_8_annot <- data.annotation.gr[data.annotation.gr$nearestGeneName %in% c("ELFN1", "LHX2", "MGMT", "PTDSS2", "ADD3", "CUEDC2", "HARS2", "PI4KB")]

# Subset raw data set by seq names. range and strand
gene_1833 <- subsetByOverlaps(gene_8, gene_8_annot, type = "equal")

# ------------ Visualization bt Heatmap ------------
gene_names <- "PI4KB"
annota_tem <- gene_8_annot[gene_8_annot$nearestGeneName == gene_names]
gene_tem <- subsetByOverlaps(gene_1833, annota_tem, type = "equal")
plot <- gene_tem %>% as.data.frame() %>% dplyr::arrange(start) %>%
  dplyr::select(-c("seqnames", "end", "width", "strand"))

tem <- plot %>% select(-"start")
plot_t <- data.frame(data.table::transpose(tem))
rownames(plot_t) <- base::colnames(plot)[-1]
base::colnames(plot_t) <- plot$start

disease <- plot_t %>%
  rownames_to_column("sample") %>%
  select(sample) %>%
  left_join(group %>% dplyr::select(SampleID, DxGroup1), 
            by = c("sample" = "SampleID"))

# ------- Export
pdf(file = "/Users/xu/OneDrive/SingleCell/Heatmap for Region/PI4KB.pdf",
    width = 15,
    height = 18)

Heatmap(plot_t, row_split = disease$DxGroup1, cluster_columns = F,
        column_title = gene_names,
        column_title_gp = gpar(fontsize = 30, fontface = "bold"),
        row_title_gp = gpar(fontsize = 20, fontface = "bold"),
        col = colorRamp2(c(0, 0.5, 1), c("blue", "#EEEEEE", "red"), space = "RGB"),
        na_col = "black",
        row_title_rot = 0,
        row_gap = unit(5, "mm"))

dev.off()

# save(list = ls(.GlobalEnv),
#      file = "/Users/xu/OneDrive/SingleCell/Heatmap for Region/Heatmap_vi.RData")











