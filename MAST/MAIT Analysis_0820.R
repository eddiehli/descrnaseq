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
})

# ------------- Pre-process
load("/Users/xu/Downloads/ipf_qc_genes_union_biogrid_intact.Robj")  # Load Seurat object into r (Name : ipf)
ipf <- subset(ipf, subset = nFeature_RNA > 950 & nFeature_RNA < 1000  & nCount_RNA < 2500) # 2626 samples
ipf <- subset(ipf, features = sample(c(1:dim(ipf[["RNA"]]@data)[1]), 500))  # Subset first 1000 features for testing : 500 features * 2626 samples
# save(ipf, file = "/Users/xu/OneDrive/SingleCell/SEURAT/Wilcox/ipd_test.Robj")

# ------------ MAST Tutorial ---------
data(maits, package = "MAST") # Load the tutorial dataset : maits (3 large lists)
dim(maits$expressionmat)
head(maits$cdat) # 96 * 21
head(maits$fdat) # 16302 * 3 gene names 

freq_expressed <- 0.2 # Set some constants
FCTHRESHOLD <- log2(1.5)  # Set some constants

scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat) # 96 single cells measured over 16302 genes
head(scaRaw)

# ---- Filtering
filterCrit <- with(colData(scaRaw), pastFastqc=="PASS"& exonRate >0.3 & PercentToHuman>0.6 & nGeneOn> 4000) # Filtering
sca <- subset(scaRaw,filterCrit)
eid <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys = mcols(sca)$entrez,keytype ="GENEID",columns = c("GENEID","TXNAME"))
ueid <- unique(na.omit(eid)$GENEID)
sca <- sca[mcols(sca)$entrez %in% ueid,]
## Remove invariant genes
sca <- sca[sample(which(freq(sca)>0), 6000),]  # 6000 features * 73 cells

# Recalculating the cellular detection rate (ngeneson)
cdr2 <- colSums(assay(sca)>0)
qplot(x=cdr2, y=colData(sca)$nGeneOn) + xlab('New CDR') + ylab('Old CDR')
colData(sca)$cngeneson <- scale(cdr2)


scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=value))+geom_density() +facet_wrap(~symbolid, scale='free_y')

thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
par(mfrow=c(5,4))
plot(thres)

# limit ourselves to genes that are found in at least 0.2 of the sample
assays(sca, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes,]  # 2349 * 73

# ----- Differential Expression using a Hurdle model
# set the reference level of the factor to be the “unstimulated” cells
cond <- factor(colData(sca)$condition)
cond <- relevel(cond, "Unstim") # relevel the unstim as the reference
colData(sca)$condition <- cond
zlmCond <- zlm(~ condition + cngeneson, sca)
slotNames(zlmCond) # Get all names of results

# log likelihood test 
lrt <- lrTest(zlmCond, "condition")
lrt <- lrTest(zlmCond, CoefficientHypothesis('conditionStim'))
dimnames(lrt) # Get names of this array
lrt[ , , 1]["1831", ]  # z=1:lambda z=2:df z=3:Pr(>Chisq)


#only test the condition coefficient.
summaryCond <- summary(zlmCond) 
#print the top 4 genes by contrast using the logFC
print(summaryCond, n = 4)
print(summaryCond, n = 4, by = "D")
print(summaryCond, n = 4, by = "C")

summaryCond$datatable %>% filter(component %in% c("C", "D")) %>% 
  filter(contrast != "(Intercept)") %>%
  arrange(desc(z))

summaryCond$datatable %>% 
  dplyr::filter(primerid == "1831") %>% 
  dplyr::filter(contrast != "(Intercept)")


summaryCond$datatable %>% 
  dplyr::filter(contrast == "conditionStim" & component == "S") %>% 
  as.data.frame() %>%
  ggplot(aes(x = z)) +
  geom_histogram(bins = 100)



