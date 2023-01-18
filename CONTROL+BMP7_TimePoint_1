# R version 4.1.3(x64 bit) and RStudio 2022.07.1+554 "Spotted Wakerobin" Release (7872775ebddc40635780ca1ed238934c3345c5de, 2022-07-22) for Windows
Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.8 Chrome/69.0.3497.128 Safari/537.36  were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/
# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed
install.packages("devtools")
install.packages("usethis")
library(usethis)
library(devtools)
install.packages("pkgbuild")
library(pkgbuild)
install.packages("Matrix")
install.packages("ggridges")
install.packages("cowplot")
install.packages('ggrepel')
install.packages("R.utils")
install.packages("gridExtra")
install.packages("Seurat")
install.packages("plotly")
install.packages("clustree")
install.packages('multtest')
install.packages("ggplot2")
install.packages("ggraph")
install.packages('Seurat')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")

# Install multtest for Seurat
BiocManager::install("multtest")
library(ggplot2)
library(cowplot)
library(Matrix)
library(ggridges)
library(ggrepel)
library(dplyr)
library(Seurat)
library(plotly)
library(clustree)
library(Seurat)
library(ggraph)
library(Matrix)
library(VGAM)
library(stats4)
library(splines)
library(DDRTree)
library(irlba)
library(BiocGenerics)
library(Biobase)

# Set global environment parameter
options(future.globals.maxSize= 89128960)
memory.limit()
memory.limit(size=21474836480)  

# Setup BMP7-ST2(Timepoint 1) Seuratobject

# OBJECT SETUP AND NORMALIZATION ####
#  Load 10X data ####
library(SeuratData)
library(SeuratObject)
library(Seurat)
BMP1st2.data <-Read10X(data.dir = "C:/Users/mad1188/Box/Sequencing Data/cellranger_count/ST2-BMP-7/filtered_feature_bc_matrix")
BMP2st2.data <-Read10X(data.dir = "C:/Users/mad1188/Box/Sequencing Data/cellranger_count/ST-2-BMP-7/filtered_feature_bc_matrix")
BMP3st2.data <-Read10X(data.dir = "C:/Users/mad1188/Box/Sequencing Data/cellranger_count/ST--2-bmp-7/filtered_feature_bc_matrix")

# Create Seurat objects ####
BMP1st2 <- CreateSeuratObject(counts = BMP1st2.data,
                              project = "BMP1st2"
)
BMP2st2 <- CreateSeuratObject(counts = BMP2st2.data,
                              project = "BMP2st2"
)
BMP3st2 <- CreateSeuratObject(counts = BMP3st2.data,
                              project = "BMP3st2"
)

head(x = colnames(x = BMP1st2))
BMP1st2 <- RenameCells(object = BMP1st2, add.cell.id = "BMP1st2")
head(x = colnames(x = BMP1st2))
head(x = colnames(x = BMP2st2))
BMP2st2 <- RenameCells(object = BMP2st2, add.cell.id = "BMP2st2")
head(x = colnames(x = BMP2st2))
head(x = colnames(x = BMP3st2))
BMP3st2 <- RenameCells(object = BMP3st2, add.cell.id = "BMP3st2")
head(x = colnames(x = BMP3st2))



# Thresholding ####
# RNA based cell thresholding
# The operator can add columns to object metadata. This is a great place to stash QC stats

BMP1st2[["percent.mt"]] <- PercentageFeatureSet(object = BMP1st2, pattern = "^MT-")
BMP2st2[["percent.mt"]] <- PercentageFeatureSet(object = BMP2st2, pattern = "^MT-")
BMP3st2[["percent.mt"]] <- PercentageFeatureSet(object = BMP3st2, pattern = "^MT-")
p1 <- VlnPlot(object = BMP1st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(object = BMP2st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(object = BMP3st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(p1, p2, p3), ncol =1)
p4 <- FeatureScatter(object = BMP1st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p5 <- FeatureScatter(object = BMP1st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p6 <- FeatureScatter(object = BMP2st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p7 <- FeatureScatter(object = BMP2st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p8 <- FeatureScatter(object = BMP3st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p9 <- FeatureScatter(object = BMP3st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p4, p5, p6, p7, p8, p9), ncol =2)
BMP1st2 <- subset(x = BMP1st2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
BMP2st2 <- subset(x = BMP2st2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
BMP3st2 <- subset(x = BMP3st2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
p10 <- VlnPlot(object = BMP1st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p11 <- VlnPlot(object = BMP2st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p12 <- VlnPlot(object = BMP3st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(p10, p11, p12), ncol =1)
p13 <- FeatureScatter(object = BMP1st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p14 <- FeatureScatter(object = BMP1st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p15 <- FeatureScatter(object = BMP2st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p16 <- FeatureScatter(object = BMP2st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p17 <- FeatureScatter(object = BMP3st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p18 <- FeatureScatter(object = BMP3st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p13, p14, p15, p16, p17, p18), ncol =2)



BMP7ST2 <- sample(x = c("BMP7ST2"), size = ncol(x = BMP1st2), replace = TRUE)
BMP1st2$groups <- BMP7ST2
head(x = BMP1st2[[]])

BMP7ST2 <- sample(x = c("BMP7ST2"), size = ncol(x = BMP2st2), replace = TRUE)
BMP2st2$groups <- BMP7ST2
head(x = BMP2st2[[]])

BMP7ST2 <- sample(x = c("BMP7ST2"), size = ncol(x = BMP3st2), replace = TRUE)
BMP3st2$groups <- BMP7ST2
head(x = BMP3st2[[]])

# Setup CONTROL-ST2(Timepoint 1) Seuratobject
# OBJECT SETUP AND NORMALIZATION ####
# Load 10X data ####
# Create Seurat objects ####

ctrl1st2.data <-Read10X(data.dir = "C:/Users/mad1188/Box/Sequencing Data/cellranger_count/ST2-CTRL/filtered_feature_bc_matrix")
ctrl1st2 <- CreateSeuratObject(counts = ctrl1st2.data,
                               project = "ctrl1st2"
)
ctrl2st2.data <-Read10X(data.dir = "C:/Users/mad1188/Box/Sequencing Data/cellranger_count/ST-2-CTRL-1/filtered_feature_bc_matrix")
ctrl2st2 <- CreateSeuratObject(counts = ctrl2st2.data,
                               project = "ctrl2st2")

ctrl3st2.data <-Read10X(data.dir = "C:/Users/mad1188/Box/Sequencing Data/cellranger_count/ST-2-CTRL-2/filtered_feature_bc_matrix")
ctrl3st2 <- CreateSeuratObject(counts = ctrl3st2.data,
                               project = "ctrl3st2"
)
head(x = colnames(x = ctrl1st2))
ctrl1st2 <- RenameCells(object = ctrl1st2, add.cell.id = "ctrl1st2")
head(x = colnames(x = ctrl1st2))
head(x = colnames(x = ctrl2st2))
ctrl2st2 <- RenameCells(object = ctrl2st2, add.cell.id = "ctrl2st2")
head(x = colnames(x = ctrl2st2))
head(x = colnames(x = ctrl3st2))
ctrl3st2 <- RenameCells(object = ctrl3st2, add.cell.id = "ctrl3st2")
head(x = colnames(x = ctrl3st2))


# Thresholding ####
# RNA based cell thresholding
# The operator can add columns to object metadata. This is a great place to stash QC stats

ctrl1st2[["percent.mt"]] <- PercentageFeatureSet(object = ctrl1st2, pattern = "^MT-")
ctrl2st2[["percent.mt"]] <- PercentageFeatureSet(object = ctrl2st2, pattern = "^MT-")
ctrl3st2[["percent.mt"]] <- PercentageFeatureSet(object = ctrl3st2, pattern = "^MT-")
p101 <- VlnPlot(object = ctrl1st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p201 <- VlnPlot(object = ctrl2st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p301 <- VlnPlot(object = ctrl3st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(p101, p201, p301), ncol =1)
p401 <- FeatureScatter(object = ctrl1st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p501 <- FeatureScatter(object = ctrl1st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p601 <- FeatureScatter(object = ctrl2st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p701 <- FeatureScatter(object = ctrl2st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p801 <- FeatureScatter(object = ctrl3st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p901 <- FeatureScatter(object = ctrl3st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p401, p501, p601, p701, p801, p901), ncol =2)
ctrl1st2 <- subset(x =  ctrl1st2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
ctrl2st2 <- subset(x = ctrl2st2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
ctrl3st2 <- subset(x = ctrl3st2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
p1001 <- VlnPlot(object = ctrl1st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1101 <- VlnPlot(object = ctrl2st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1201 <- VlnPlot(object = ctrl3st2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(p1001, p1101, p1201), ncol =1)
p1301 <- FeatureScatter(object = ctrl1st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1401 <- FeatureScatter(object = ctrl1st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1501 <- FeatureScatter(object = ctrl2st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1601 <- FeatureScatter(object = ctrl2st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1701 <- FeatureScatter(object = ctrl3st2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1801 <- FeatureScatter(object = ctrl3st2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p1301, p1401, p1501, p1601, p1701, p1801), ncol =2)



ctrlST2 <- sample(x = c("ctrlST2"), size = ncol(x = ctrl1st2), replace = TRUE)
ctrl1st2$groups <- ctrlST2
head(x = ctrl1st2[[]])

ctrlST2 <- sample(x = c("ctrlST2"), size = ncol(x = ctrl2st2), replace = TRUE)
ctrl2st2$groups <- ctrlST2
head(x = ctrl2st2[[]])

ctrlST2 <- sample(x = c("ctrlST2"), size = ncol(x = ctrl3st2), replace = TRUE)
ctrl3st2$groups <- ctrlST2
head(x = ctrl3st2[[]])


# Merge Datasets
# Based on comment to Issue #4753 https://github.com/satijalab/seurat/issues/4753
# We use RPCA to yield conserved mapping

DiabetesST2<- c(ctrl1st2, ctrl2st2, ctrl3st2,BMP1st2,BMP2st2,BMP3st2)

# normalize and identify variable features for each dataset independently
DiabetesST2.list <- lapply(X = DiabetesST2, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = DiabetesST2.list)
DiabetesST2.list <- lapply(X = DiabetesST2.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Perform integration
DiabetesST2.anchors <- FindIntegrationAnchors(object.list = DiabetesST2.list, anchor.features = features, reduction = "rpca")
DiabetesST2.combined <- IntegrateData(anchorset = DiabetesST2.anchors)

# Run the standard workflow for visualization and clustering
DefaultAssay(DiabetesST2.combined) <- "integrated"
DiabetesST2.combined <- ScaleData(DiabetesST2.combined, verbose = FALSE)
DiabetesST2.combined <- RunPCA(DiabetesST2.combined, npcs = 30, verbose = FALSE)
DiabetesST2.combined <- RunUMAP(DiabetesST2.combined, reduction = "pca", dims = 1:30)
DiabetesST2.combined <- FindNeighbors(DiabetesST2.combined, reduction = "pca", dims = 1:30)
DiabetesST2.combined <- FindClusters(DiabetesST2.combined, resolution = 0.4)

# Visualization
p3 <- DimPlot(DiabetesST2.combined, reduction = "umap", raster = FALSE)
p4 <- DimPlot(DiabetesST2.combined, reduction = "umap", label = TRUE,raster = FALSE,
              repel = TRUE)
p3 + p4

# Linear dimensionality assessment
DefaultAssay(object = DiabetesST2.combined)
DefaultAssay(object = DiabetesST2.combined) <- "integrated"
DiabetesST2.combined <- RunPCA(object = DiabetesST2.combined, features = VariableFeatures(object = DiabetesST2.combined))
print(x = DiabetesST2.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = DiabetesST2.combined, dims = 1:2, reduction = "pca")
DimPlot(object = DiabetesST2.combined, reduction = "pca")
DimHeatmap(object = DiabetesST2.combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = DiabetesST2.combined, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(object = DiabetesST2.combined)

#CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
DiabetesST2.combined <- FindNeighbors(object = DiabetesST2.combined, dims = 1:20)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.1)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.2)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.3)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.4)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.5)
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(DiabetesST2.combined, prefix = "integrated_snn_res.")
DiabetesST2.combined <- FindClusters(object = DiabetesST2.combined, resolution = 0.4)

table(DiabetesST2.combined$seurat_clusters)
table(Idents(DiabetesST2.combined), DiabetesST2.combined$orig.ident)
DiabetesST2.combined <- RunUMAP(DiabetesST2.combined, dims = 1:30)
p1010 <- DimPlot(object = DiabetesST2.combined, group.by = c("orig.ident", "integrated_snn_res.0.4"), combine = FALSE, pt.size = 1, raster = FALSE, label = TRUE)
p1010 <- lapply(X = p1010, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(p1010)

# DATA-VISUALIZATION ####
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DiabetesST2UMAP.combined <- RenameIdents(DiabetesST2.combined, `0` = "A2.1", `1` = "D1", `2` = "A5", `3` = "D3", `4` = "A4", `5` = "A2", `6` = "D6", `7` = "A7", `8` = "A1", `9` = "S",`10` = "D2", `11` = "M", `12` = "DX", `13` = "I", `14` = "D4", `15` = "E", `16` = "D7")
DimPlot(DiabetesST2UMAP.combined, label = TRUE, raster = FALSE)
DimPlot(DiabetesST2UMAP.combined, label = FALSE, raster = FALSE)

DimPlot(DiabetesST2UMAP.combined, reduction = "umap",
        cols =  c("sienna",
                  "turquoise1","yellow",
                  "royalblue1",
                  "greenyellow",
                  "sienna3",
                  "darkturquoise",
                  "gold3",
                  "lightgoldenrod2",
                  "grey",
                  "mediumseagreen",
                  "pink",
                  "darkcyan",
                  "red",
                  "lightslateblue",
                  "palegreen2",
                  "navyblue"
        ),
        pt.size = 1, raster = FALSE, label = TRUE)

#Renaming the clusters
DiabetesST2UMAP.combined_AD <- RenameIdents(DiabetesST2.combined, `1` = "D1",`10` = "D2",`3` = "D3",`14` = "D4",`6` = "D6",`16` = "D7",`12` = "DX",`8` = "A1",`5` = "A2",`0` = "A2.1", `4` = "A4", `2` = "A5",`7` = "A7", `9` = "S", `11` = "M",  `13` = "I",  `15` = "E")
DimPlot(DiabetesST2UMAP.combined_AD, label = TRUE, raster = FALSE)
DimPlot(DiabetesST2UMAP.combined_AD, label = FALSE, raster = FALSE)

DimPlot(DiabetesST2UMAP.combined_AD, reduction = "umap",
        cols =  c("turquoise1","mediumseagreen",
                  "royalblue1","lightslateblue",
                  "darkturquoise","navyblue", "darkcyan",
                  "lightgoldenrod2","sienna3",
                  "sienna","greenyellow",
                  "yellow","gold3",
                  "grey",
                  "pink",
                  "red",
                  "palegreen2"
        ),
        pt.size = 1, raster = FALSE, label = TRUE)

DiabetesST2UMAP.combined_AD$CellType <- Idents(DiabetesST2UMAP.combined_AD)

DefaultAssay(object = DiabetesST2UMAP.combined_AD) <- "RNA"
DiabetesST2UMAP.combined_AD <- NormalizeData(DiabetesST2UMAP.combined_AD)
DiabetesST2UMAP.combined_AD <- SCTransform(DiabetesST2UMAP.combined_AD, assay = "RNA", new.assay.name = "SCT", verbose = TRUE, return.only.var.genes = TRUE)
DefaultAssay(object = DiabetesST2.combined)
DefaultAssay(object = DiabetesST2.combined) <- "SCT"
top10.genes.vst <- head(x = VariableFeatures(object = DiabetesST2.combined), 10)
p106 <- VariableFeaturePlot(object = DiabetesST2.combined, assay = 'SCT', selection.method = c('sct'))
p106
LabelPoints(plot = p106, points = top10.genes.vst, repel = TRUE)


head(x = colnames(x = DiabetesST2.combined))
tail(x = colnames(x = DiabetesST2.combined))

unique(x = sapply(X = strsplit(x = colnames(x = DiabetesST2.combined), split = "_"), FUN = "[", 1))
table(DiabetesST2.combined$orig.ident)

DefaultAssay(object = DiabetesST2.combined)
DefaultAssay(object = DiabetesST2.combined) <- "SCT"
Diabetes9.markers <- FindAllMarkers(object = DiabetesST2.combined,
                                    features = VariableFeatures(DiabetesST2.combined, assay = 'SCT'),
                                    only.pos = TRUE,
                                    min.pct = 0.1,
                                    logfc.threshold = 0.41,
                                    assay = 'SCT',
                                    slot = c('data'))
library(magrittr)
library(dplyr)

#Finding differentially expressed features (cluster biomarkers)
library(magrittr)
library(dplyr)
DiabetesST2.combined.markers <- FindAllMarkers(DiabetesST2.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DiabetesST2.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(DiabetesST2.combined.markers, 'C:/Users/mad1188/Box/Diabetes_new/The whole UMAP/BMPmarkers.csv')

# save seuratobject
saveRDS(DiabetesST2UMAP.combined, file = "C:/Users/mad1188/Box/Mayur_shared_JDB/JUAN ONLY-DO NOT TOUCH/Mayur last/Materials and Methods/CONTROL+BMP7_TimePoint_1.rds")
