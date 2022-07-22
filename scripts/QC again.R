# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "RawData/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "azimuth_pbmc", min.cells = 10, min.features = 200) 
#min.cells: Include features detected in at least this many cells.

# Compute percent mito ratio
pbmc$mitoRatio <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmc$mitoRatio<-pbmc@meta.data$mitoRatio/100

# Create metadata dataframe
metadata_pbmc <- pbmc@meta.data

# Add cell IDs to metadata
metadata_pbmc$cells <- rownames(metadata_pbmc)

# Rename columns
library(tidyverse)
metadata_pbmc <- metadata_pbmc %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
pbmc@meta.data <- metadata_pbmc

#filter at cell level-------------------------------------------------------------------------------------
# Filter out low quality reads using selected thresholds 
max_nGene<-mean(pbmc@meta.data$nGene)+4*sd(pbmc@meta.data$nGene) #5765.344
filtered_pbmc <- subset(x = pbmc, 
                        subset= (nGene <= max_nGene) &
                          (nGene >= 200) & 
                          (mitoRatio <= 0.1))

#filter at gene level-------------------------------------------------------------------------------------
# Extract counts
counts <- GetAssayData(object = filtered_pbmc, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_PBMC<- CreateSeuratObject(filtered_counts, meta.data = filtered_pbmc@meta.data)

# Normalize the counts
filtered_PBMC1 <- NormalizeData(filtered_PBMC)

# Score cells for cell cycle
filtered_PBMC1 <- CellCycleScoring(filtered_PBMC1, 
                                   g2m.features = cc.genes$g2m.genes, 
                                   s.features = cc.genes$s.genes)

# run sctransform
filtered_PBMC1 <- SCTransform(filtered_PBMC1, vars.to.regress = c("mitoRatio"), verbose = F)


# These are now standard steps in the Seurat workflow for visualization and clustering
filtered_PBMC1 <- RunPCA(filtered_PBMC1, npcs=30,verbose = F)
ElbowPlot(filtered_PBMC1)
filtered_PBMC1 <- FindNeighbors(filtered_PBMC1, dims = 1:30, verbose = F)
filtered_PBMC1 <- FindClusters(filtered_PBMC1, verbose = F)
filtered_PBMC1 <- RunUMAP(filtered_PBMC1, dims = 1:30, verbose = F)
filtered_PBMC1 <- RunTSNE(filtered_PBMC1, dims = 1:30, verbose = F)
DimPlot(filtered_PBMC1,
        reduction = "umap",
        label = TRUE,
        label.size = 5)

DimPlot(filtered_PBMC1,
        reduction = "tsne",
        label = TRUE,
        label.size = 5)
