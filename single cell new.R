# -------------------------------
# 📦 Load Required Libraries
# -------------------------------
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SingleR", quietly = TRUE)) BiocManager::install("SingleR")
if (!requireNamespace("celldex", quietly = TRUE)) BiocManager::install("celldex")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")

library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
library(SummarizedExperiment)

# -------------------------------
# 📥 Load 10x Genomics Data
# -------------------------------
data_dir <- "D:/filtered_feature_bc_matrix/"
sc_data <- Read10X(data.dir = data_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "scRNAseq_Project", min.cells = 3, min.features = 200)

# -------------------------------
# 🧼 Quality Control
# -------------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visual QC
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# -------------------------------
# 🧪 Normalization and HVGs
# -------------------------------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Visualize top 10 HVGs
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# -------------------------------
# ⚖️ Scaling and PCA
# -------------------------------
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Elbow plot for PC selection
ElbowPlot(seurat_obj)

# Estimate expected number of doublets
nExp <- round(0.075 * ncol(seurat_obj))  # 7.5% is a common estimate

# Run DoubletFinder
seurat_obj <- doubletFinder_v3(seurat_obj,
                               PCs = 1:10,
                               pN = 0.25,
                               pK = 0.09,  # Needs to be optimized
                               nExp = nExp,
                               reuse.pANN = FALSE,
                               sct = FALSE)

# -------------------------------
# 🧱 Clustering and UMAP
# -------------------------------
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

# -------------------------------
# 🧬 Find Cluster Markers
# -------------------------------
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top_markers$gene) + NoLegend()

# -------------------------------
# 🧠 Automated Cell Type Annotation with SingleR
# -------------------------------
ref <- celldex::HumanPrimaryCellAtlasData()
norm_data <- GetAssayData(seurat_obj, slot = "data")
singleR_results <- SingleR(test = norm_data, ref = ref, labels = ref$label.main)

# Add SingleR labels to metadata
seurat_obj$SingleR.labels <- singleR_results$labels

# Plot UMAP with SingleR annotations
DimPlot(seurat_obj, group.by = "SingleR.labels", label = TRUE, reduction = "umap") +
  ggtitle("SingleR Cell Type Annotation")

# -------------------------------
# 🔬 Differential Gene Expression (DGE) Analysis
# -------------------------------

# Set identities to clusters and compare cluster 0 vs 1
Idents(seurat_obj) <- "seurat_clusters"
dge_clusters <- FindMarkers(seurat_obj, ident.1 = 0, ident.2 = 1, logfc.threshold = 0.25, min.pct = 0.25)
write.csv(dge_clusters, "D:/scRNA_project/DGE_cluster0_vs_cluster1.csv")

# Set identities to SingleR annotations and compare B cells vs T cells
Idents(seurat_obj) <- seurat_obj$SingleR.labels
dge_celltypes <- FindMarkers(seurat_obj, ident.1 = "B cells", ident.2 = "T cells", logfc.threshold = 0.25)
write.csv(dge_celltypes, "D:/scRNA_project/DGE_B_vs_T_cells.csv")

# -------------------------------
# 💾 Save Results
# -------------------------------
saveRDS(seurat_obj, file = "D:/scRNA_project/seurat_annotated.rds")
write.csv(cluster_markers, "D:/scRNA_project/cluster_markers.csv")
