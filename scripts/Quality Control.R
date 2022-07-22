#install the latest version of Seurat
install.packages('Seurat')
packageVersion("Seurat") #‘4.1.0’; Depends: R (≥ 4.0.0)；my R version 4.0.5 
#QC standard was according to Cytoimmgen paper


setwd("~/Summer2022/Azimuth_PBMC")

library(Seurat)
library(ggplot2)
library(sctransform)
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

pdf("fig/Vln_nGene&nUMImitoRatio.pdf", width = 14, height = 12)
VlnPlot(pbmc, features = c("nGene", "nUMI", "mitoRatio"), ncol = 3,pt.size =0.05)
dev.off()


pdf("fig/FeatureScatter_QC.pdf", width = 10, height = 8)
FeatureScatter(object = pbmc, feature1 = "nUMI", feature2 = "mitoRatio")
FeatureScatter(object = pbmc, feature1 = "nUMI", feature2 = "nGene")
dev.off()

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


##reassess QC matrix-------------------------------------------------------------------------------------
## Visualize QC metrics as a violin plot
pdf("fig/Vln_nGene&nUMImitoRatio_afterQC.pdf", width = 14, height = 12)
VlnPlot(filtered_PBMC, features = c("nGene", "nUMI", "mitoRatio"), ncol = 3,pt.size =0.05)
dev.off()

metadata_clean <- filtered_PBMC@meta.data

pdf("fig/nGene&nCell_afterQC.pdf", width = 10, height = 7)
# Visualize the number of cell counts per sample
p1<-metadata_clean %>% 
  ggplot(aes(x=seq_folder, fill=seq_folder)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells")

# Visualize the distribution of genes detected per cell via boxplot
p2<-metadata_clean %>% 
  ggplot(aes(x=seq_folder, y=nGene, fill=seq_folder)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nGene per cell")

p1+p2
dev.off()

pdf("fig/density_afterQC.pdf", width = 16, height = 5)
# Visualize the number UMIs/transcripts per cell
p3<-metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=nUMI, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 

# Visualize the distribution of genes detected per cell via histogram
p4<-metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=nGene, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(200,max_nGene),linetype="dashed")

# Visualize the distribution of mitochondrial gene expression detected per cell
p5<-metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=mitoRatio, fill=seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.1, linetype="dashed")
p3+p4+p5
dev.off()

pdf("fig/correlation_afterQC.pdf", width = 10, height = 7)
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_hline(yintercept = c(200,max_nGene), linetype="dashed")
dev.off()


# Normalize the counts
filtered_PBMC1 <- NormalizeData(filtered_PBMC)

# Score cells for cell cycle
filtered_PBMC1 <- CellCycleScoring(filtered_PBMC1, 
                                   g2m.features = cc.genes$g2m.genes, 
                                   s.features = cc.genes$s.genes)

# run sctransform
filtered_PBMC1 <- SCTransform(filtered_PBMC1, vars.to.regress = c("mitoRatio","S.Score", "G2M.Score"), verbose = F)


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


#doublet deletion---------------------------------------------------------------------------------------
library(DoubletFinder)


# pK Identification (no ground-truth) 
sweep.res.list_PBMC <- paramSweep_v3(filtered_PBMC1, PCs = 1:30, sct = T)
sweep.stats_PBMC <- summarizeSweep(sweep.res.list_PBMC, GT = FALSE)
bcmvn_PBMC <- find.pK(sweep.stats_PBMC)

ggplot(bcmvn_PBMC,aes(pK,BCmetric,group=1))+
  geom_point()+
  geom_line()+
  theme(axis.text.x=element_text(angle = 45))

pK<-bcmvn_PBMC %>%
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  select(pK)

pK<-as.numeric(as.character(pK[[1]]))

# Homotypic Doublet Proportion Estimate 
annotations <- filtered_PBMC1@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.076*nrow(filtered_PBMC1@meta.data))  #7.6% is from https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies 
filtered_PBMC1 <- doubletFinder_v3(filtered_PBMC1, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
#pN:of artificial doublets, default value 0.25 but authors have found DoubletFinder performance is not dependent on this value
#pK: the neighborhood size (pK) used to compute the number of artificial nearest neighbors, DoubletFinder performance dependent on this value

#visualize doublets
names(filtered_PBMC1@meta.data)
DimPlot(filtered_PBMC1,reduction = "umap",group.by = "DF.classifications_0.25_0.05_564")
VlnPlot(filtered_PBMC1, features = "nGene", group.by = "DF.classifications_0.25_0.05_564", pt.size = 0.1)

#number of singlet and doublet---------------------------------------------------------------------------------------
table(filtered_PBMC1@meta.data$DF.classifications_0.25_0.05_564)

metadata <- filtered_PBMC1@meta.data
names(metadata)
colnames(metadata)[17] <- "doublet_finder"
filtered_PBMC1@meta.data <- metadata 

# subset and save
filtered_PBMC.singlets <- subset(filtered_PBMC1, doublet_finder == "Singlet")##whether do I need to rescale de subset data?
filtered_PBMC1
filtered_PBMC.singlets
DimPlot(filtered_PBMC.singlets,reduction = "umap")

saveRDS(filtered_PBMC.singlets, "metadata/filtered_PBMC_singlets_withoutregression.rds")

# featureplot the second density peak
metadata_clean <- filtered_PBMC.singlets@meta.data
metadata_clean %>% 
  ggplot(aes(color=seq_folder, x=nGene, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(200,2580,max_nGene),linetype="dashed")

FeaturePlot(pbmc, 
            reduction = "umap", 
            features ="G2M.Score", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = F)

pbmc4<-pbmc1
pbmc4<-subset(pbmc4,subset=(nGene>2580))
pdf("fig/pbmc_UMAP_secondpeak.pdf",width = 10,height =7)
DimPlot(pbmc4, reduction = "umap",  cols = c(brewer.pal(5, "Reds")[1:5],brewer.pal(4, "Blues")[1:4],brewer.pal(7, "Greens")[2:3],brewer.pal(7, "Greens")[6:7], brewer.pal(3, "Purples")[2:3], brewer.pal(9, "Oranges")[c(4,6,8)],brewer.pal(8, "Dark2")[c(2,4,8)]),
        label = TRUE,label.size = 2,repel = TRUE)
dev.off()


# vlnplot of cell type and phase number
VlnPlot(pbmc1, features = c("mitoRatio"),pt.size = 0)





