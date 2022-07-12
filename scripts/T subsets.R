setwd("~/Summer2022/Azimuth_PBMC")
library(Seurat)
library(ggplot2)

#subset T cells(have already subsetted Ts)
pbmc_T <- FindNeighbors(object = pbmc_T, 
                        dims = 1:40)
pbmc_T <- FindClusters(object = pbmc_T,
                       resolution = 3)
pbmc_T <- RunUMAP(pbmc_T, dims = 1:40, verbose = F)
DimPlot(pbmc_T,
        reduction = "umap",
        label = TRUE,
        label.size = 3)

#find more CD4s markers
FeaturePlot(pbmc_T, 
            reduction = "umap", 
            features ="AC004585.1", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = F)
#TRDC, TRGC1, TRGC2, KLRC1, NKG7, TRDV2, CD7, TRGV9, KLRD1, KLRG1

pbmc_T <- FindNeighbors(object = pbmc_T, 
                      dims = 1:40)
pbmc_T <- FindClusters(object = pbmc_T,
                     resolution = 4.5)
Idents(object = pbmc_T) <- "SCT_snn_res.4.5"
DimPlot(pbmc_T,label=T)

saveRDS(pbmc_T,"metadata/pbmc_T.rds")

pbmc_T1<-pbmc_T
# Rename all identities
pbmc_T1 <- RenameIdents(object = pbmc_T1, 
                      "0" = "CD8 Naive",
                      "1" = "CD4 Naive",
                      "2" = "CD4 TCM",
                      "3" = "CD4 TEM",
                      "4" = "CD4 TCM",
                      "5" = "MAIT",
                      "6" = "CD4 Naive",
                      "7" = "CD4 Naive",
                      "8" = "CD4 Naive",
                      "9" = "CD4 TCM",
                      "10" = "Treg",
                      "11" = "CD4 Naive",
                      "12" = "CD4 Naive",
                      "13" = "CD4 TCM",
                      "14" = "CD4 TCM",
                      "15" = "CD4 TCM",
                      "16" = "CD4 TCM",
                      "17" = "CD4 Naive",
                      "18" = "CD4 Naive",
                      "19" = "CD8 TCM",
                      "20" = "CD4 TCM",
                      "21" = "CD4 Naive",
                      "22" = "CD4 TCM",
                      "23" = "CD8 TEFF",
                      "24" = "CD8 TEM",
                      "25" = "CD4 TCM",
                      "26" = "CD4 TEFF",
                      "27" = "CD4 TCM",
                      "28" = "CD4 TCM",
                      "29" = "CD4 TCM",
                      "30" = "gdT",
                      "31" = "CD4 Naive")
DimPlot(pbmc_T1,label=T)
Idents(pbmc_T1)<-factor(Idents(pbmc_T1),levels = c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEFF","MAIT","gdT"))
library(RColorBrewer)
colour_T1 = c(brewer.pal(5, "Reds")[1:5],brewer.pal(5, "Blues")[1:5],brewer.pal(7, "Greens")[3])
pdf("fig/T_UMAP_labled.pdf",width = 10,height =7)
DimPlot(pbmc_T1,label=T,cols = colour_T1,pt.size = 0.5,label.size = 3,repel = T)
dev.off()

CD4Naive<-WhichCells(pbmc_T1,idents = "CD4 Naive")
CD4TCM<-WhichCells(pbmc_T1,idents = "CD4 TCM")
CD4TEM<-WhichCells(pbmc_T1,idents = "CD4 TEM")
Treg<-WhichCells(pbmc_T1,idents = "Treg")
CD4TEFF<-WhichCells(pbmc_T1,idents = "CD4 TEFF")
CD8Naive<-WhichCells(pbmc_T1,idents = "CD8 Naive")
CD8TCM<-WhichCells(pbmc_T1,idents = "CD8 TCM")
CD8TEM<-WhichCells(pbmc_T1,idents = "CD8 TEM")
CD8TEFF<-WhichCells(pbmc_T1,idents = "CD8 TEFF")
MAIT<-WhichCells(pbmc_T1,idents = "MAIT")
gdT<-WhichCells(pbmc_T1,idents = "gdT")
pbmc2<-pbmc1
metatable<-pbmc2@meta.data[,c(8,24)]
metatable$active.ident<-factor(metatable$active.ident,
                               levels = c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEFF","MAIT","gdT","XCL1 NK","NK","NK Proliferating","Naive B","Memory B","CD14 Mono","CD16 Mono","mDC","pDC","MGK","HSPC"))
metatable$active.ident[which(metatable$cells %in% CD4Naive)]<-"CD4 Naive"
metatable$active.ident[which(metatable$cells %in% CD4TCM)]<-"CD4 TCM"
metatable$active.ident[which(metatable$cells %in% CD4TEM)]<-"CD4 TEM"
metatable$active.ident[which(metatable$cells %in% Treg)]<-"Treg"
metatable$active.ident[which(metatable$cells %in% CD4TEFF)]<-"CD4 TEFF"
metatable$active.ident[which(metatable$cells %in% CD8Naive)]<-"CD8 Naive"
metatable$active.ident[which(metatable$cells %in% CD8TCM)]<-"CD8 TCM"
metatable$active.ident[which(metatable$cells %in% CD8TEM)]<-"CD8 TEM"
metatable$active.ident[which(metatable$cells %in% CD8TEFF)]<-"CD8 TEFF"
metatable$active.ident[which(metatable$cells %in% MAIT)]<-"MAIT"
metatable$active.ident[which(metatable$cells %in% gdT)]<-"gdT"

pbmc2@meta.data[,c(8,24)]<-metatable
colour_all = c(brewer.pal(5, "Reds")[1:5],brewer.pal(5, "Blues")[1:5],brewer.pal(7, "Greens")[2:3],brewer.pal(7, "Greens")[6:7], brewer.pal(3, "Purples")[2:3], brewer.pal(9, "Oranges")[c(4,6,8)],brewer.pal(8, "Dark2")[c(2,4,8)])
pdf("fig/pbmc_UMAP_newlabled.pdf",width = 10,height =7)
DimPlot(pbmc2,label=T,group.by = "active.ident",cols = colour_all,repel = T,label.size = 2)
     +theme(plot.title = element_blank())
dev.off()
