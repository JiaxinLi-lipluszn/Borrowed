# Process the perturb seq

# Build Seurat object from expression matrix

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/jiaxinli/Desktop/data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

mydata <- read_csv("/Users/jiaxinli/downloads/perturb_10_classes.csv")
cell_name <- mydata$X1
mydata <- mydata[,-1]
row.names(mydata) <- cell_name
trans <- t(mydata)
mydataMatrix <- as(as.matrix(trans), "dgCMatrix")

mydata_S <- CreateSeuratObject(counts = mydataMatrix, project = "perturb", min.cells = 3, min.features = 200)
pbmc

mydata_S <- NormalizeData(mydata_S) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()



DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')



bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

# export scaled data for RNA (2000 variable genes) and protein
rna.scale <- GetAssayData(bm[["RNA"]], slot="scale.data")[VariableFeatures(bm[["RNA"]]), ]
protein.scale <- GetAssayData(bm[["ADT"]], slot="scale.data")
metadata <- bm@meta.data

print(dim(rna.scale))
print(dim(protein.scale))

write.csv(rna.scale, gzfile("data/rna_scale.csv.gz"))
write.csv(protein.scale, gzfile("data/protein_scale.csv.gz"))
write.csv(metadata, gzfile("data/metadata.csv.gz"))
