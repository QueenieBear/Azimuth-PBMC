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
            features ="NCR3", 
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




#compare with Azimuth prediction -----------------------------------------------------------------------------
pbmc_T1$active.ident<-factor(pbmc_T1@active.ident,levels = c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEFF","MAIT","gdT"))
saveRDS(pbmc_T1,"metadata/pbmc_T1.rds")

#pbmc_T1 <- readRDS("~/Summer2022/Azimuth_PBMC/metadata/pbmc_T1.rds")
predictions <- read.delim('metadata/azimuth_pred.tsv', row.names = 1)
pbmc_T1 <- AddMetaData(
  object =pbmc_T1 ,
  metadata = predictions)

colnames(pbmc_T1@meta.data)
predictable<-table(pbmc_T1@meta.data$active.ident,pbmc_T1@meta.data$predicted.celltype.l2)
library(ztable)
pdf("fig/heatmapTable.pdf",width = 6.5,height = 5)
pheatmap(predictable, 
         scale = "none", 
         cluster_cols = F, 
         cluster_rows = F,
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         main="Compared with Azimuth Prediction(x-axis)",
         border_color = NA,
         display_numbers = T,
         number_format = "%.0f",
         number_color="white",
         legend=F,
         legend_breaks=c(0,10,20,50,100,200,400,1000,2000),
         color = c(brewer.pal(9, "Blues")[9:5]),
         breaks =c(0,10,20,50,100,200,400,1000,2000),
         cellwidth=20,
         cellheight=20)
dev.off()





metatable2<-pbmc_T1@meta.data[,c(8,24,25)]

metatable2$mismatched[which(metatable2$active.ident != metatable2$predicted.celltype.l2)] <-1
metatable2$mismatched[which(metatable2$active.ident == metatable2$predicted.celltype.l2)] <-0
metatable2<-metatable2[which(metatable2$mismatched == "1"),]
metatable2$RealMismatched[which((metatable2$active.ident=="CD4 Naive")&(metatable2$predicted.celltype.l2=="CD4 TCM"))]<-1
metatable2$RealMismatched[which((metatable2$active.ident=="CD4 TEM")&(metatable2$predicted.celltype.l2=="CD4 TCM"))]<-1
metatable2$RealMismatched[which((metatable2$active.ident=="CD4 TEM")&(metatable2$predicted.celltype.l2=="CD8 TCM"))]<-1
metatable2$RealMismatched[which((metatable2$active.ident=="Treg")&(metatable2$predicted.celltype.l2=="CD4 TCM"))]<-1
metatable2$RealMismatched[which((metatable2$active.ident=="CD8 TCM")&(metatable2$predicted.celltype.l2=="CD8 TEM"))]<-1

metatable2<-metatable2[which(metatable2$RealMismatched == "1"),]
metatable2$newname<-paste("mismatched",metatable2$active.ident)
metatable2$newname<-paste(metatable2$newname,"/")
metatable2$newname<-paste(metatable2$newname,metatable2$predicted.celltype.l2)
metatable2<-metatable2[,c(1,6)]
pbmc_T1@meta.data<-left_join(pbmc_T1@meta.data,metatable2,by="cells")
#colnames(pbmc_T1@meta.data)
metatable3<-pbmc_T1@meta.data[,c(8,24,28)]
metatable3$active.ident<-as.character(metatable3$active.ident)
metatable3$newcelltype<-ifelse(is.na(metatable3$newname),metatable3$active.ident,metatable3$newname)
pbmc_T1@meta.data<-left_join(pbmc_T1@meta.data,metatable3[c(1,4)],by="cells")
a<-table(pbmc_T1@meta.data$newcelltype)
b<-as.factor(names(a))
c<-levels(b)
pbmc_T1@meta.data$newcelltype<-factor(pbmc_T1@meta.data$newcelltype,levels =c)
pbmc_T1$newcelltype

mismatchFeature<-c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","TRAC","TRBC2","TRBC1",
                   "CCR7","FHIT","LEF1","RGS10","TCF7","NOSIP","SELL","IL7R","PASK","PEACAM1",
                   "IL7R","TMSB10","ITGB1","LTB","AQP3","MAL","CCL5","GZMK","GZMH","FGFBP2","GNLY","NKG7","LINC02446","ANXA1","CXCR4","KLRB1","GZMA","TXN","TNFRSF4","TNFRSF1A","TNFRSF1B","TNFRSF6B","	TNFRSF8","LYAR",
                   "NKG7","GNLY","CCL4","CCL5","GZMB","GZMH","IL32",
                   "CD44","CD69",
                   "RTKN2","FOXP3","AC133644.2","IL2RA","CCL5","KLRD1","GZMK","CST7","TRGC2","TIGIT","CTLA4",
                   "KLRB1","NKG7","GZMK","SLC4A10","NCR3",
                   "TRDC","TRGC1","KLRC1","NKG7")
mismatchFeature<-unique(mismatchFeature)
pdf("fig/mismatch_dotplot.pdf", width = 18, height = 7)
DotPlot(pbmc_T1, features = mismatchFeature, dot.scale = 6,group.by = "newcelltype") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 4, name = "RdBu"))) + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,family = "Times",size = rel(0.7)))
dev.off()


pbmc_T2<-pbmc_T
metatable4<-pbmc_T1@meta.data[,c(8,29)]
metatable4$newcelltype[which(metatable4$newcelltype=="mismatched CD4 Naive / CD4 TCM")]<-"CD4 Naive"
metatable4$newcelltype[which(metatable4$newcelltype=="mismatched CD4 TEM / CD4 TCM")]<-"CD4 Naive"
metatable4$newcelltype[which(metatable4$newcelltype=="mismatched CD4 TEM / CD8 TCM")]<-"CD4 TEM"
metatable4$newcelltype[which(metatable4$newcelltype=="mismatched CD8 TCM / CD8 TEM")]<-"CD8 TEM"
metatable4$newcelltype[which(metatable4$newcelltype=="mismatched Treg / CD4 TCM")]<-"Treg"
metatable4$newcelltype<-factor(metatable4$newcelltype,levels = c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 TEFF","CD8 Naive","CD8 TCM","CD8 TEM","CD8 TEFF","MAIT","gdT"))

pbmc_T2@meta.data[,c(8,25)]<-metatable4

pdf("fig/UMAP_T_unlabeled.pdf",width = 10,height = 7)
DimPlot(pbmc_T2,label=F,cols = colour_T1,group.by = "newcelltype",pt.size = 0.5,label.size = 3,repel = T)+ggtitle("Newly Defined Celltypes")
DimPlot(pbmc_T1,label=F,cols = colour_T1,group.by = "active.ident",pt.size = 0.5,label.size = 3,repel = T)+ggtitle("Previously Defined Celltypes")
dev.off()


pdf("fig/UMAP_T_labeled.pdf",width = 10,height = 7)
DimPlot(pbmc_T2,label=T,cols = colour_T1,group.by = "newcelltype",pt.size = 0.5,label.size = 3,repel = T)+ggtitle("Newly Defined Celltypes")
DimPlot(pbmc_T1,label=T,cols = colour_T1,group.by = "active.ident",pt.size = 0.5,label.size = 3,repel = T)+ggtitle("Previously Defined Celltypes")
dev.off()













