setwd("./")
source('libraryAndFunctions.R')
dir.create('figure')
projectName<-'TEC'
CPM.data<-read.csv('./CPMForSeurat.csv',
                   row.names = 1)

## Preprocessing
mca<-preprocessData(data=CPM.data, 
                    metadata = NULL,
                    projectName = projectName,
                    normalization.method = 'LogNormalize', # normalization of normal methods
                    scale.factor = 1000, # normalization paramater of normal methods
                    selection.methods = "vst",
                    nfeatures=2000,
                    doTsneTransform = TRUE,  # whether do t-SNE transformation
                    perp = 20,  # perplexity
                    fin_dims = 30, # final dimension of transform space
                    ini_dims = NULL,
                    platform = 'MATLAB',
                    slot = 'data',
                    selectedFeatures = 300
)

## Visulization with TSNE
mca<-RunTSNE(mca,reduction = 'tsnetransform',dims = 1:30,reduction.name = 'tsne_tsnet')

## KNN-Clustering
mca<-FindNeighbors(mca,dims = 1:30, reduction = 'tsnetransform',k.param = 20,graph.name = 'graph_tsnet')
mca<-FindClusters(mca,graph = 'graph_tsnet',resolution = 0.38)

levels(mca)<-c('6','4','3','0','1','2','7','5')
mca<-RenameIdents(mca,
                  '6' = 'Cluster 1',
                  '4' = 'Cluster 2',
                  '3' = 'Cluster 3',
                  '0' = 'Cluster 4',
                  '1' = 'Cluster 5',
                  '2' = 'Cluster 6',
                  '7' = 'Cluster 7',
                  '5' = 'Cluster 8')

##Scatter Plot
h<-DimPlot(mca,reduction='tsne_tsnet',pt.size=3)
h<-h+theme(legend.text = element_text(size = 30),legend.key.size = unit(1.5,'cm'),
           legend.title=element_blank(),
           axis.title = element_text(size = 20))
ggsave(filename = './figure/ScatterCluster.png',
       plot = h,
       height = 20,
       width = 24,
       units = 'cm',
       dpi = 600,
       limitsize = FALSE,
       scale = 1.2
)

h<-DimPlot(mca,reduction='tsne_tsnet',group.by = 'orig.ident',cols = c('orangered1','goldenrod'),pt.size=3)
h<-h+theme(legend.text = element_text(size = 30),legend.key.size = unit(1.5,'cm'),
           legend.title=element_blank(),
           axis.title = element_text(size = 20))
ggsave(filename = './figure/DayCluster.png',
       plot = h,
       height = 20,
       width = 24,
       units = 'cm',
       dpi = 600,
       limitsize = FALSE,
       scale = 1.2
)


for(i in c('Foxn1','Psmb11','Notch1','Notch2','Notch3','Dll4','Jag1','Jag2','Mfng','Hes6','Heyl','Foxn1','Psmb11')){
  h<-plot_density(mca,feature=i,reduction='tsne_tsnet',size=3)
  ggsave(filename = paste0('./densityScatter_',i,'.png'),
         plot = h,
         height = 20,
         width = 24,
         units = 'cm',
         dpi = 600,
         limitsize = FALSE,
         scale = 1.2
  )
}

##Plot Heatmap
markers<- FindAllMarkers(object = mca, only.pos = TRUE, min.pct = 0, assay = 'RNA',
                         logfc.threshold = 0.25)
markers %>%  group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC) -> top
h<-DoHeatmap(mca, features = top$gene ,assay = 'RNA')
h<-h+scale_fill_gradient2(high = 'brown2',mid ='lightcyan',low = 'steelblue2')
ggsave(filename = './figure/ClusterHeatmap.png',
       plot = h,
       height = 15,
       width = 12,
       units = 'cm',
       dpi = 1200,
       limitsize = FALSE,
       scale = 2
)

##VlnPlot
targetmarkers<-c('Mt1','Nefm','Tagln3','Neurod1','Epcam','Fgf8','Plet1','Pcna','Pax1','Foxn1','Prss16','Ccl25','Krt18','Krt8','Cldn3','Cldn4','Krt5','Krt17','Krt19','Pth','Id4')
myColor<-c('#f8766d','#cd9600','#7cae00','#00be67','#00bfc4','#00a9ff','#c77cff','#ff61cc')
h<-multiVlnPlot(mca,markerGene=targetmarkers,myColorCluster = myColor)
ggsave(filename = './figure/ClusterVln.png',
       plot = h,
       height = 10,
       width = 35,
       units = 'cm',
       dpi = 1200,
       limitsize = FALSE,
       scale = 1
)


##===================================================
##Regulon analysis with SCENIC
library(SCENIC)
scenicOptions <- initializeScenic(org="mgi", dbDir="../../H/SCENIC/", nCores=2)
exprMat<-GetAssayData(mca,slot = 'data')
exprMat<-as.matrix(exprMat)
#exprMat<-exprMat[deg,]
runCorrelation(exprMat, scenicOptions)
runGenie3(exprMat, scenicOptions)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
#scenicOptions@settings$nCores=2
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat)
save.image('outputkm.RData')
cellInfo <- data.frame(seuratCluster=Idents(mca))
runSCENIC_4_aucell_binarize(scenicOptions)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType <- regulonActivity_byCellType[rowMax(regulonActivity_byCellType)>0.025,]
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
save.image('outputkm.RData')

png(filename = '/media/sf_Projects/figure/Thymus/Poster/heatmapSCENIC_loga.png',
    height = 50,
    width = 25,
    units = 'cm',
    res = 1200)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity (%)", col = c("chartreuse4","wheat","violetred4")
                        ,row_names_max_width = ComplexHeatmap::max_text_width(rownames(regulonActivity_byCellType_Scaled),gp=gpar(fontsize=20,fontface='bold'))
                        ,width = unit(10, "cm")
                        ,column_names_gp = gpar(fontsize = 20,fontface="bold")
                        ,row_names_gp = gpar(fontsize = 20,fontface="bold")
                        ,heatmap_legend_param = list(legend_height = unit(100,'mm')
                                                     , legend_weight = unit(80,'mm')
                                                     , title_position = 'leftcenter-rot'
                                                     , title_gp = gpar(fontsize = 30, fontface = "bold")
                                                     , labels_gp = gpar(fontsize = 20)
                                                     
                        )
)
dev.off()

regulons<-loadInt(scenicOptions,'regulons')
regulons('Sox9')


library(SCENIC)
dir.create('./SCENIC')
setwd('./SCENIC')
data(list="motifAnnotations_mgi_v9",package="RcisTarget")
motifAnnotations_mgi<-motifAnnotations_mgi_v9

# Download mm9-500bp-upstream-7species.mc9nr.feather and mm9-tss-centered-10kb-7species.mc9nr.feather from resouces.aertslab.org
scenicOptions <- initializeScenic(org="mgi", dbDir="./", nCores=8)
exprMat<-GetAssayData(mca,slot = 'data')
exprMat<-as.matrix(exprMat)
runCorrelation(exprMat, scenicOptions)
runGenie3(exprMat, scenicOptions)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat)
cellInfo <- data.frame(seuratCluster=Idents(mca))
runSCENIC_4_aucell_binarize(scenicOptions)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType <- regulonActivity_byCellType[rowMax(regulonActivity_byCellType)>0.02,]
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
library(ComplexHeatmap)
png(filename = '../figure/heatmapSCENIC.png',
    height = 50,
    width = 25,
    units = 'cm',
    res = 1200)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity (%)", col = c("chartreuse4","wheat","violetred4")
                        ,row_names_max_width = ComplexHeatmap::max_text_width(rownames(regulonActivity_byCellType_Scaled),gp=gpar(fontsize=20,fontface='bold'))
                        ,width = unit(10, "cm")
                        ,column_names_gp = gpar(fontsize = 20,fontface="bold")
                        ,row_names_gp = gpar(fontsize = 20,fontface="bold")
                        ,heatmap_legend_param = list(legend_height = unit(100,'mm')
                                                     , legend_weight = unit(80,'mm')
                                                     , title_position = 'leftcenter-rot'
                                                     , title_gp = gpar(fontsize = 30, fontface = "bold")
                                                     , labels_gp = gpar(fontsize = 20)
                                                     
                        )
)
dev.off()
