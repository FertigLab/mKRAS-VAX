library(Seurat)
library(parallel)
library(dplyr)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
setwd("~/salvage_tmp/final/")
set.seed(99)
load(file = './code_outputs/2024-05-18_AggrSeurat.rda')
batch1 = subset(seu, subset = Batch == 1)
batch2 = subset(seu, subset = Batch == 2)
rm(seu)

#Process Batch1============================================================
pbmc = batch1; rm(batch1)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
plot0 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- plot1 + plot2
plot3

filt <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 7.5)
ncells_orig = ncol(pbmc); rm(pbmc)

plot4 <- FeatureScatter(filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot6 <- plot4 + plot5
plot6
plot7 <- VlnPlot(filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells_remaining = ncol(filt)

filt <- NormalizeData(filt, normalization.method = "LogNormalize", scale.factor = 10000)
filt <- FindVariableFeatures(filt, selection.method = 'vst', nfeatures = 5000)

top10 <- head(VariableFeatures(filt), 10)
plot8 <- VariableFeaturePlot(filt)
plot9 <- LabelPoints(plot = plot8, points = top10, repel = TRUE)
plot9

filt <- ScaleData(filt, features = VariableFeatures(filt))
filt <- RunPCA(filt, features = VariableFeatures(filt))
plot10 <- ElbowPlot(filt, ndims = 30)

N=25
filt<- RunUMAP(filt, dims = 1:N)
filt = filt[sample(c(1:nrow(filt)), nrow(filt), replace = F),] #shuffle cells
plot11 <- DimPlot(filt, reduction = 'umap', group.by = 'Sample', raster = FALSE)
plot12 <- DimPlot(filt, reduction = 'umap',  raster = FALSE, split.by = 'Sample')




umap = filt@reductions$umap@cell.embeddings
filt@meta.data$to_exclude = umap[,1] < -6 & umap[,2] < 0
filt2 = subset(filt, to_exclude == FALSE)
plot13 <- DimPlot(filt2, reduction = 'umap', group.by = 'Sample', raster = FALSE)
plot14 <- DimPlot(filt2, reduction = 'umap',  raster = FALSE, split.by = 'Sample')



filt2 <- FindNeighbors(filt2, dims = 1:25) #low nDim because I'm just trying to catch the batch-specific cluster 
filt2 <- FindClusters(filt2, resolution = 0.3)
save(filt2, file = paste0('./code_outputs/', Sys.Date(), '_batch1filt.rda'))

#Process Batch2============================================================
pbmc = batch2; rm(batch2)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
plot0 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- plot1 + plot2
plot3

filt <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 7.5)
ncells_orig = ncol(pbmc); rm(pbmc)

plot4 <- FeatureScatter(filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot6 <- plot4 + plot5
plot6
plot7 <- VlnPlot(filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ncells_remaining = ncol(filt)

filt <- NormalizeData(filt, normalization.method = "LogNormalize", scale.factor = 10000)
filt <- FindVariableFeatures(filt, selection.method = 'vst', nfeatures = 5000)

top10 <- head(VariableFeatures(filt), 10)
plot8 <- VariableFeaturePlot(filt)
plot9 <- LabelPoints(plot = plot8, points = top10, repel = TRUE)
plot9

filt <- ScaleData(filt, features = VariableFeatures(filt))
filt <- RunPCA(filt, features = VariableFeatures(filt))
plot10 <- ElbowPlot(filt, ndims = 50)

N=40
filt<- RunUMAP(filt, dims = 1:N)
filt = filt[sample(c(1:nrow(filt)), nrow(filt), replace = F),] #shuffle cells
plot11 <- DimPlot(filt, reduction = 'umap', group.by = 'Sample', raster = FALSE)
plot12 <- DimPlot(filt, reduction = 'umap',  raster = FALSE, split.by = 'Sample')

#Based on UMAP let's try patient correction using Batchelor. 
library(batchelor)
counts <- filt@assays$RNA@counts
counts <- counts[rownames(counts) %in% VariableFeatures(filt),]
counts = as.matrix(counts)
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
sce = SingleCellExperiment(assays = list(logcounts = log2(t(t(counts)/size.factors) + 1)))
meta = filt@meta.data
#reducedDims(sce) <- list(PCA = filt@reductions$pca@cell.embeddings, UMAP = filt@)
batch = filt@meta.data$Patient == 'J1994.002'
rm(counts)
set.seed(99)
sce = fastMNN(
  sce, 
  batch = batch
)

#Jerryrig solution to use native seurat functions on Batchelor dimensional reduction
pca = filt@reductions$pca@cell.embeddings
tmp = reducedDim(sce, "corrected")
colnames(tmp) = colnames(pca)
filt@reductions$corrected = filt@reductions$pca
filt@reductions$corrected@cell.embeddings = tmp

filt <- RunUMAP(filt, reduction.name = 'umap_corrected', reduction = 'corrected', dims = 1:N)

plot13 <- DimPlot(filt, reduction = 'umap_corrected', group.by = 'Sample', raster = FALSE)
plot14 <- DimPlot(filt, reduction = 'umap_corrected',  raster = FALSE, split.by = 'Sample')


filt <- FindNeighbors(filt, dims = 1:40) 
filt <- FindClusters(filt, resolution = 0.3)
save(filt, file = paste0('./code_outputs/', Sys.Date(), '_batch2filt.rda'))