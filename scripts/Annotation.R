setwd("~/Summer2022/Azimuth_PBMC")

library(Seurat)
pbmc <- readRDS("~/Summer2022/Azimuth_PBMC/metadata/filtered_PBMC_singlets.rds")

#marker from https://www.sciencedirect.com/science/article/pii/S1074761320303162?via%3Dihub#mmc1
#CD3G, KLRF1, XCL1,MS4A1,IGHG1,MZB1,CD68,LYZ,MKI67,TOP2A,GZMA,PPBP

intermB<-c("TNFRSF13B", "IGHM", "IGHD","AIM2", "LINC01857", "RALGPS2")
memB<-c("COCH", "AIM2", "BANK1", "SSPN",  "TEX9", "RALGPS2", "TNFRSF13C", "LINC01781")
naiveB<-c("IGHM", "IGHD",  "IL4R", "MS4A1", "CXCR4", "BTG1", "TCL1A",  "YBX3")
plasmablast<-c("IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5", "HRASLS2", "NT5DC2")
platelet<-c("PPBP", "PF4", "NRGN", "GNG11", "CAVIN2", "TUBB1", "CLU", "HIST1H2AC", "RGS18", "GP9")
NKprolif<-c("MKI67", "KLRF1", "TYMS", "TRDC", "TOP2A", "FCER1G", "PCLAF", "CD247", "CLSPN", "ASPM")
CD8TCM<-c("ANXA1", "KRT1", "LINC02446", "YBX3", "IL7R", "TRAC", "NELL2", "LDHB")
CD8TEM<-c("CCL5", "GZMH", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "TRGC2")
HSPC<-c("SPINK2", "PRSS57", "CYTL1", "EGFL7", "GATA2", "CD34", "SMIM24", "AVP", "MYB", "LAPTM4B")
NK_CD56bright<-c("XCL2", "FCER1G", "SPINK2", "TRDC", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1", "TNFRSF11A")
HIVpaper<-c("CD3D","CD4","CCR7","SELL","TCF7","LEF1","S100A4","GPR183","TNFRSF4",
              "GZMK","GZMA","GZMB","GZMH","GNLY","NKG7",
              "KLRD1","TRDC","TYROBP","FCGR3A","CD38",
              "TIGIT","HAVCR2","CTLA4","LAG3","PDCD1",
              "IFI6","ISG15","MKI67","PCNA","SLC4A10","TRAV1-2",
              "TRGV9","TRDV2","FOXP3","IL2RA")
#annotation by marker genes
FeaturePlot(pbmc, 
            reduction = "umap", 
            features ="CD8B", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = F)



# clustering---------------------------------------------------------------------------------------------

pbmc <- RunPCA(pbmc, npcs=40,verbose = F)

                              

DefaultAssay(pbmc)<-"SCT"
# Determine the clusters for various resolutions  
pbmc <- FindNeighbors(object = pbmc, 
                      dims = 1:40)
pbmc <- FindClusters(object = pbmc,
                     resolution = 4.5)

# Plot the UMAP
DimPlot(pbmc,
        reduction = "umap",
        label = TRUE,
        label.size = 4)

pbmc1<-pbmc
# Rename all identities
pbmc1 <- RenameIdents(object = pbmc1, 
                     "0" = "CD4 TEM",
                     "1" = "NK",
                     "2" = "CD14 Mono",
                     "3" = "CD14 Mono",
                     "4" = "CD14 Mono",
                     "5" = "CD14 Mono",
                     "6" = "CD14 Mono",
                     "7" = "CD4 Naive",
                     "8" = "Naive B",
                     "9" = "CD4 TEM",
                     "10" = "CD8 Naive",
                     "11" = "CD4 Naive",
                     "12" = "CD4 Naive",
                     "13" = "CD14 Mono",
                     "14" = "CD14 Mono",
                     "15" = "Naive B",
                     "16" = "CD14 Mono",
                     "17" = "CD4 TCM",
                     "18" = "CD4 TCM",
                     "19" = "CD4 Naive",
                     "20" = "CD4 Naive",
                     "21" = "CD16 Mono",
                     "22" = "CD4 Naive",
                     "23" = "CD4 TCM",
                     "24" = "CD14 Mono",
                     "25" = "MAIT",
                     "26" = "CD8 TEM",
                     "27" = "CD8 TEM",
                     "28" = "CD14 Mono",
                     "29" = "Treg",
                     "30" = "Memory B",
                     "31" = "mDC",
                     "32" = "CD14 Mono",
                     "33" = "Memory B",
                     "34" = "XCL1 NK",
                     "35" = "Naive B",
                     "36" = "CD8 TEFF",
                     "37" = "Memory B",
                     "38" = "NKT",
                     "39" = "MGK",
                     "40" = "CD4 TEFF",
                     "41" = "HSPC",
                     "42" = "NK Proliferating",
                     "43" = "pDC")


library(RColorBrewer)
library(ggforce)
Idents(pbmc1)<-factor(Idents(pbmc1),levels = c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TEM","CD8 TEFF","MAIT","NKT","XCL1 NK","NK","NK Proliferating","Naive B","Memory B","CD14 Mono","CD16 Mono","mDC","pDC","MGK","HSPC"))

pdf("fig/pbmc_UMAP_labled.pdf",width = 10,height =7)
DimPlot(pbmc1, reduction = "umap",  cols = c(brewer.pal(5, "Reds")[1:5],brewer.pal(4, "Blues")[1:4],brewer.pal(7, "Greens")[2:3],brewer.pal(7, "Greens")[6:7], brewer.pal(3, "Purples")[2:3], brewer.pal(9, "Oranges")[c(4,6,8)],brewer.pal(8, "Dark2")[c(2,4,8)]),
        label = TRUE,label.size = 2,repel = TRUE)
dev.off()


saveRDS(pbmc,"metadata/pbmc.rds")
saveRDS(pbmc1,"metadata/pbmc1.rds")



#marker genes display---------------------------------------------------------------------------------------------

pbmc_T<-subset(pbmc1,idents = c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TEM","CD8 TEFF","MAIT","NKT"))
pbmc2<-pbmc1
pbmc2 <- RenameIdents(object = pbmc2, 
                      "CD4 Naive"= "T",
                      "CD4 TCM"= "T",
                      "CD4 TEM"= "T",
                      "Treg"= "T",
                      "CD4 TEFF"= "T",
                      "CD8 Naive"= "T",
                      "CD8 TEM"= "T",
                      "CD8 TEFF"= "T",
                      "MAIT"= "T",
                      "NKT"= "T")
colour_2 = c(brewer.pal(4, "Blues")[4],brewer.pal(7, "Greens")[3],brewer.pal(7, "Greens")[6:7], brewer.pal(3, "Purples")[2:3], brewer.pal(9, "Oranges")[c(4,6,8)],brewer.pal(8, "Dark2")[c(2,4,8)])
Markers_2<-c("CD3D","XCL1","KLRF1","KLRC1","KLRD1","TRDC","MKI67","PCNA","TOP2A","PCLAF","CLSPN",
             "MS4A1","CD79A","CD19","TCL1A","IGHD","IL4R","TCL1A","IGHG1","COCH","AIM2","SSPN",
             "LYZ","CD68","CD14","FCGR3A","CD1C","CLEC4C","LILRA4","IL3RA",
             "PPBP","CMTM5","ITGB3","MYL9",
             "SPINK2", "PRSS57", "CYTL1", "EGFL7", "GATA2")

source("scripts/StackedVlnPlot.R")
pdf("fig/vlnplot_Markers_2.pdf", width = 5, height = 25)
StackedVlnPlot(pbmc2, features = Markers_2, pt.size = 0, cols = colour_2)
dev.off()

colour_T = c(brewer.pal(5, "Reds")[1:5],brewer.pal(4, "Blues")[1:4],brewer.pal(7, "Greens")[1])
Markers_T<-c("CD3D","CCR7","SELL","LEF1","S100A4","GPR183","TNFRSF4","CD4","FOXP3","IL2RA","CD8A","CD8B","GZMH","GNLY","NKG7","GZMK","GZMA","GZMB","SLC4A10","NCR3")

source("scripts/StackedVlnPlot.R")
pdf("fig/vlnplot_Markers_T.pdf", width = 4, height = 14)
StackedVlnPlot(pbmc_T, features = Markers_T, pt.size = 0, cols = colour_T)
dev.off()



#top10 genes display---------------------------------------------------------------------------------------------
DefaultAssay(pbmc1)<-"SCT"

# Connect to AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()

head(unique(ah$species))
unique(ah$species)

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)



# Find markers for every cluster compared to all remaining cells, report only the positive ones
allmarkers<-FindAllMarkers(object = pbmc1, 
                           only.pos = TRUE,
                           logfc.threshold = 0.25) 

# Extract top 10 markers per cluster
top10 <- allmarkers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, 
        wt = avg_log2FC)

top10_anno<-top10 %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name")) %>%
  cbind(cluster_id = top10$cluster, .)

# View data
write.csv(top10_anno,file="metadata/top10_anno.csv")

#heatmap
pdf("fig/heatmap_top10.pdf", width = 18, height = 15)
colour_all = c(brewer.pal(5, "Reds")[1:5],brewer.pal(4, "Blues")[1:4],brewer.pal(7, "Greens")[2:3],brewer.pal(7, "Greens")[6:7], brewer.pal(3, "Purples")[2:3], brewer.pal(9, "Oranges")[c(4,6,8)],brewer.pal(8, "Dark2")[c(2,4,8)])
DoHeatmap(pbmc1, features = top10$gene,size = 2,group.colors=colour_all,draw.lines=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))+
  theme(axis.text.y = element_text(hjust = 1, vjust = 1,family = "Times",size = rel(0.7)))
dev.off()


#upload in Azimuth Model
levels(pbmc1$SCT_snn_res.4.5)
pbmc1$active.ident<-factor(pbmc1@active.ident,levels=c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TEM","CD8 TEFF","MAIT","NKT","XCL1 NK","NK","NK Proliferating","Naive B","Memory B","CD14 Mono","CD16 Mono","mDC","pDC","MGK","HSPC"))
saveRDS(pbmc1,"metadata/pbmc1.rds")

AzimuthResults <- read_csv("metadata/AzimuthResults.csv")

