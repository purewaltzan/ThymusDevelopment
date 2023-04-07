##Load Required Packages
#library('rstudioapi')
library('Seurat')
library('ggplot2')
library('dplyr')
library('pROC')
library('tsne')
library('matlabr')
library('monocle')

##=========================
## Set default tool paths
tsnetlibpath<-'.'
scentlibpath<-'.'

##=========================
## tSNE transformation 
# Class for tSNE transformation results
#dim.reduction <- setClass(
#  Class = "DimReduc",
#  slots = c(
#    cell.embeddings = 'matrix',
#    feature.loadings = 'matrix',
#    feature.loadings.projected = 'matrix',
#    assay.used = 'character',
#    stdev = 'numeric',
#    key = 'character',
#    jackstraw = 'JackStrawData',
#    misc = 'list'
#  )
#)

# tSNE transformation function
tsneTransform<-function(object,
                        perp=30, # perplexity
                        fin_dims=30, # final dimension
                        ini_dims=30, # initial dimension
                        platform='R', # the platform running t-SNE, MATLAB is recommanded
                        reduction = "tsnetransform",
                        assay = NULL,
                        weight.by.var = TRUE,
                        verbose = TRUE,
                        reduction.key = "tsnet_",
                        seed.use = 42,
                        approx = TRUE,
                        outputfile = NULL,
                        deleteTemp = FALSE,
                        slot = 'data',
                        ...){
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  cells <- colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  raw_data<-as.matrix(GetAssayData(object,assay = 'tsnet'))
  if(length(object@assays$tsnet@var.features)>10){
    raw_data<-raw_data[object@assays$tsnet@var.features,]
  }
  if(platform=='outputfiles'){
    if(is.null(outputfile)){
      outputfile<-paste(tsnetlibpath,'/tsnetoutput.csv',sep='')
      tsner<-read.csv(outputfile,header = FALSE, row.names = NULL)
    }
  }else{
    if(platform=='R'){
      warning("Trying tSNE transformation in R. If there are more than 2000 cells, it will take very long time. Try to use platform='MATLAB'")
      tsner<-tsne(t(raw_data),
                  perplexity = perp,
                  initial_dims = ini_dims,
                  k = fin_dims,
      )
    }else{
      if(platform=='MATLAB'){
        get_matlab()
        #tsnetlibpath<-'../../tsneTransformation'
        if(have_matlab()){
          write.csv(t(raw_data),file=paste(tsnetlibpath,'/data.csv',sep=''))
          write.csv(c(perp,fin_dims,ini_dims),file=paste(tsnetlibpath,'/para.csv',sep=''))
          run_matlab_script(fname=paste(tsnetlibpath,'/tsnet.m',sep=''), verbose = TRUE, desktop = FALSE,
                            splash = FALSE, display = FALSE, wait = TRUE,
                            single_thread = FALSE)
          outputfile<-paste(tsnetlibpath,'/tsnetoutput.csv',sep='')
          tsner<-read.csv(outputfile,header = FALSE, row.names = NULL)
          #deleteTemp = TRUE
        }else{
 
          warning('There is no MATLAB detected in this environment, attempting to use output files')
          if(is.null(outputfile)){
            outputfile<-paste(tsnetlibpath,'/tsnetoutput.csv',sep='')
          }
          tsner<-try(read.csv(outputfile,header = FALSE, row.names = NULL))
          if(!is.array(tsner)){
            stop("There is no t-SNE transformation output files. Please locate it or generate a new one or use platform = 'R'")
          }else{
            deleteTemp = FALSE
          }
        }
        if(deleteTemp){
          file.remove(paste(tsnetlibpath,'/data.csv',sep=''))
          file.remove(paste(tsnetlibpath,'/para.csv',sep=''))
          file.remove(paste(tsnetlibpath,'/tsnetoutput.csv',sep=''))
        }
      }else{
        stop('t-SNE transform should be run in R or MATLAB')
      }
    }
  }
  cna<-vector()
  for(i in 1:ncol(tsner)){
    cna<-c(cna,paste('tsnet_',i,sep=''))
  }
  print('1')
  colnames(tsner)<-cna
  print('2')
  rownames(tsner)<-cells
  ttf.obj <- CreateDimReducObject(embeddings = as.matrix(tsner),assay='tsnet')
#  ttf.obj <- new(Class = "DimReduc",cell.embeddings = as.matrix(tsner))
  object@reductions[[reduction]]<-ttf.obj
 # DefaultAssay(object@reductions[[reduction]])<-"tsnet"
  return(object)
}

# tSNE impletation of Seurat
preprocessData<-function(data, # expression matrix
                         metadata = NULL,
                         projectName = 'Seurat',## parameter of normal Seurat methods
                         normalization.method = 'LogNormalize', # normalization of normal methods
                         scale.factor = 1000, # normalization paramater of normal methods
                         selection.methods = "vst",
                         nfeatures=2000,## parameter of t-SNE transformation
                         doTsneTransform = FALSE,  # whether do t-SNE transformation
                         perp = 30,  # perplexity
                         fin_dims = 30, # final dimension of transform space
                         ini_dims = NULL, # dimension before transform space
                         platform = 'MATLAB',
                         slot = 'data',
                         selectedFeatures = (ncol(data)*0.9),
                         n.neighbors = 30L,...
){
  mca <- CreateSeuratObject(counts = data, 
                            meta.data = metadata,
                            project = projectName)
  if(is.null(ini_dims)){
    ini_dims<-fin_dims
  }
  if(doTsneTransform){
    n<-rowSums(data>0)
    mca@assays$tsnet<-mca@assays$RNA
    mca@assays$tsnet@var.features<-rownames(data)[n<(selectedFeatures)]
    mca <- tsneTransform(object = mca, 
                         perp = perp, 
                         fin_dims = fin_dims, 
                         ini_dims = ini_dims,
                         platform = platform,
                         slot = slot)
    mca <- RunUMAP(object = mca, 
                   dims=1:fin_dims,
                   reduction='tsnetransform',
                   reduction.name = 'umap_tsnet',reduction.key = 'UMAP_',n.neighbors = n.neighbors)
  }
  mca <- NormalizeData(object = mca, 
                       normalization.method = normalization.method, 
                       scale.factor = scale.factor)
  mca <- FindVariableFeatures(mca,
                              selection.method = "vst",nfeatures = nfeatures)
  mca <- ScaleData(object = mca)
  mca <- RunPCA(object = mca, 
                pc.genes = mca@var.genes)
  mca <- JackStraw(object = mca, 
                   num.replicate = 100)
  mca <- RunUMAP(object = mca, 
                 dims=1:20,n.neighbors = n.neighbors)
  mca <- ScaleData(object = mca,features = rownames(mca))
  return(mca)
}


GSEAroc<-function(mca,proc = NA){
  seq<-levels(mca$seurat_clusters)
  cluster<-mca$seurat_clusters
  data<-GetAssayData(mca)
  r_roc<-data.frame(row.names = (rownames(mca)))
  if(!is.na(proc)){
    library("foreach")
    library("doParallel")
    cl<- makeCluster(min(proc,detectCores()-1)) 
    registerDoParallel(cl)
    for(i in 1:length(seq)){
      r_roc[i]<-paraROC(data=data,resp=(cluster==seq[i]))
    }
    stopCluster(cl)
    }else{
      for(i in 1:length(seq)){
        r_roc<-cbind(r_roc,nonParaROC(data=data,resp=(cluster==seq[i])))
      }
    }
  colnames(r_roc)<-seq
  return(r_roc)
  }

paraROC<-function(data,resp){
    trewq<-vector()
    trewq <- foreach(j = 1:nrow(data)) %dopar% 
      pROC::auc(pROC::roc(predictor = data[j,],response = resp,direction = '<',quiet = TRUE))
    return(trewq)
}

nonParaROC<-function(data,resp){
  trewq<-vector()
  for(j in 1:nrow(data)){
    trewq[j]<-auc(roc(predictor = data[j,],
                      response = resp,
                      direction = '<',
                      quiet = TRUE
    ))
    print(j)
  }

  return(trewq)
}

enumGSEA<-function(rocs,
                   ref,
                   ont = 'BP',
                   keyType = 'SYMBOL',
                   pAdjustMethod = 'none',
                   minGSSize=10,
                   maxGSSize = 100000,
                   pvalueCutoff = Inf ){
  geoList<-list()
  for(i in 1:ncol(rocs)){
    geneList <- rocs[,i]
    names(geneList) <- rownames(rocs)
    geneList <- sort(geneList,decreasing = TRUE)
    geoList <- c(geoList,gseGO(geneList = geneList,
                 ont = ont,
                 OrgDb = ref,
                 keyType = keyType,
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize,
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = pAdjustMethod))
  }
  names(geoList)<-colnames(rocs)
  return(geoList)
}

enumGO<-function(markers,
                 mca,
                 ref,
                 ont = 'BP',
                 keyType = 'SYMBOL',
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1,
                 minGSSize=10,
                 maxGSSize = 100000,
                 readable = TRUE){
  goList<-list()
  coln<-levels(markers$cluster)
  for(i in 1:length(coln)){
    gene <- markers[markers$cluster==coln[i],'gene']
    goList[i] <- enrichGO(gene, 
                                OrgDb = ref, 
                                keyType = keyType, 
                                ont = "BP",
                                pvalueCutoff = pvalueCutoff, 
                                pAdjustMethod = pAdjustMethod, 
                                universe = rownames(mca),
                                qvalueCutoff = qvalueCutoff, 
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize,
                                readable = FALSE)
  }
  names(goList)<-levels(Idents(mca))
  return(goList)
}

getGSEAData<-function(geoList,slots='NES',na.value=0){
  GSEAData <- data.frame()
  GSEAData <- merge(GSEAData,geoList[[1]]@result[,c('ID','Description',slots)],all=TRUE,all.x=TRUE,all.y=TRUE)
  for(i in 1:length(geoList)){
    d <- geoList[[i]]@result[,c('ID','Description',slots)]
    GSEAData <- merge(GSEAData,d,all=TRUE,all.x=TRUE,all.y=TRUE)
    colnames(GSEAData)[i+2]<-names(geoList)[i]
  }
  GSEAData[is.na(GSEAData)]<-na.value
  rownames(GSEAData) <- GSEAData$ID
  return(GSEAData)
}

plotGSEAHeatmap<-function(ES){
  GOName <- factor(rep(ES$Description,ncol(ES)-2),levels=rev(ES$Description))
  groupAnnotation <- rep(colnames(ES)[-(1:2)],rep(nrow(ES),ncol(ES)-2))
  orde <- rep(nrow(ES):1,ncol(ES)-2)
  data <- ES[,-(1:2)]
  print(dim(data))
  ESl <- vector()
  for(i in 1:ncol(data)){
    for(j in 1:nrow(data))
      ESl <- rbind(ESl,data[j,i])
  }
  print(length(ESl))
  print(length(groupAnnotation))
  kkdata <-data.frame(orde,ES=ESl,GOName,groupAnnotation)

  h <- ggplot(data = kkdata, aes(x = groupAnnotation, y = GOName)) +
    geom_tile(aes(fill = ES)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90)
    )
  return(h)
}

multiVlnPlot<-function(mca,markerGene,myColorCluster,vlnAdjust=2){
  gene <- factor(rep(markerGene,ncol(mca)),levels=markerGene)
  cell <- factor(rep(colnames(mca),rep(length(markerGene),ncol(mca))),levels=colnames(mca))
  cluster <- factor(rep(Idents(mca),rep(length(markerGene),ncol(mca))),levels=rev(levels(Idents(mca))))
  data <- GetAssayData(mca)[markerGene,]
  d<-matrix(data,ncol=1)
  d<-d+rnorm(n=length(d))/100000
  kkdata <- data.frame(expression = d, gene = gene, cell = cell, cluster = cluster)
  h <- ggplot(data = kkdata, aes(y = expression, x = cluster))+
    geom_violin(scale='width',aes(fill = cluster),trim = F,adjust=vlnAdjust) + 
    facet_wrap(~gene,ncol=length(markerGene),nrow=1,scales = "free_x") + 
    scale_fill_manual(values=rev(myColorCluster))+
    theme_bw()+
    coord_flip()+
    theme(legend.position = 'none',panel.border = element_rect(size = 2),
          
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(face = 'bold',size=15),
          axis.text.y = element_text(face = 'bold',size=15,colour = 'black'))
  
}


plotGSEAGODot<-function(ES,p){
  GOName <- factor(rep(ES$Description,ncol(ES)-2),levels=rev(ES$Description))
  groupAnnotation <- factor(rep(colnames(ES)[-(1:2)],rep(nrow(ES),ncol(ES)-2)),levels=colnames(ES))
  orde <- rep(nrow(ES):1,ncol(ES)-2)
  data <- ES[,-(1:2)]
  p.data<-p[,-(1:2)]
  ESl <- vector()
  pl <- vector()
  for(i in 1:ncol(data)){
    for(j in 1:nrow(data)){
      ESl <- rbind(ESl,data[j,i])
      pl <- rbind(pl,p.data[j,i])
    }
  }
  pl=-log10(pl)
  print(length(ESl))
  print(length(groupAnnotation))
  kkdata <-data.frame(orde,NES=ESl,logp=pl,GOName,groupAnnotation)
  
  h <- ggplot(data = kkdata, aes(x = groupAnnotation, y = GOName, size=(logp+0.1), colour=NES)) +
    geom_point() +
    scale_size(trans = 'log10') +
    scale_colour_gradient2('NES',high = 'brown2',mid = 'wheat',low = 'royalblue4') +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90),
      axis.text= element_text(size=10,color = 'black')
    )
  h$labels$size<-'-log(p)'
  h$labels$colour<-'NES'
  
  return(h)
}

ploCnet<-function(geneList,filepre, OrgDb, showCategory = 10){
  geneList = sort(geneList, decreasing = TRUE)
  ego <- gseGO(geneList = geneList[geneList!=0],
               OrgDb        = OrgDb,
               keyType      = "SYMBOL",
               ont          = "BP",
               nPerm        = 100,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               verbose      = FALSE)
  ego@result<-ego@result[order(ego@result$enrichmentScore,decreasing = TRUE),]
  h<-cnetplot(ego, foldChange=geneList,showCategory = showCategory)
  ggsave(filename=paste("./figure/",filepre,"up_cnet.png",sep=''),
         plot = h,
         width = 30,
         height = 30,
         unit = 'cm',
         dpi=1000,
         scale = 0.9,
         limitsize = FALSE)
  geneList = -geneList
  geneList = sort(geneList, decreasing = TRUE)
  ego <- gseGO(geneList     = geneList[geneList!=0],
               OrgDb        = OrgDb,
               keyType      = "SYMBOL",
               ont          = "BP",
               nPerm        = 100,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               verbose      = FALSE)
  ego@result<-ego@result[order(ego@result$enrichmentScore,decreasing = TRUE),]
  h<-cnetplot(ego, foldChange=-geneList,showCategory = showCategory)
  ggsave(filename=paste("./figure/",filepre,"down_cnet.png",sep=''),
         plot = h,
         width = 30,
         height = 30,
         unit = 'cm',
         dpi=1000,
         scale = 0.9,
         limitsize = FALSE)
}
ploDot<-function(geneList,filepre,OrgDb,showCategory = 20){
  geneList = sort(geneList, decreasing = TRUE)
  ego <- gseGO(geneList     = geneList[geneList!=0],
               OrgDb        = OrgDb,
               keyType      = "SYMBOL",
               ont          = "BP",
               nPerm        = 100,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               verbose      = FALSE)
  ego@result<-ego@result[order(ego@result$enrichmentScore,decreasing = TRUE),]
  h<-dotplot(ego ,showCategory = showCategory)
  ggsave(filename=paste("./figure/",filepre,"up_dot.png",sep=''),
         plot = h,
         width = 50,
         height = 20,
         unit = 'cm',
         dpi=1000,
         scale = 0.8,
         limitsize = FALSE)
  geneList=-geneList
  geneList = sort(geneList, decreasing = TRUE)
  ego <- gseGO(geneList     = geneList[geneList!=0],
               OrgDb        = OrgDb,
               keyType      = "SYMBOL",
               ont          = "BP",
               nPerm        = 100,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               verbose      = FALSE)
  ego@result<-ego@result[order(ego@result$enrichmentScore,decreasing = FALSE),]
  h<-dotplot(ego ,showCategory = showCategory)
  ggsave(filename=paste("./figure/",filepre,"down_dot.png",sep=''),
         plot = h,
         width = 50,
         height = 20,
         unit = 'cm',
         dpi=1000,
         scale = 0.8,
         limitsize = FALSE)
}

scent<-function(datak,libpath=scentlibpath){
  load('./code/ppiAsigH-PC2-17Jan2016.Rd')
  maskr<-intersect(rownames(datak),rownames(ppiAsigHMC.m))
  write.csv(x = datak[maskr,], file = paste(libpath,'data.csv',sep = ''))
  write.csv(x = ppiAsigHMC.m[maskr,maskr],file = paste(libpath,'PPI.csv',sep = ''))
  remove(ppiAsigHclass.m)
  remove(ppiAsigHMC.m)
  remove(ppiAsigHclass.v)
  remove(ppiAsigHMC.v)
}

tsnetROC<-function(mca_1,mca_2,gene){
  roc_log_auc<-0
  roc_tsnet_auc<-0
  for(i in levels(mca_1)){
    temp_roc<-roc(Idents(mca_1)==i,GetAssayData(mca_1,'counts')[gene,],quiet = TRUE)
    if((temp_roc$auc>roc_log_auc)&(sum(Idents(mca_1)==i)>5)){
      roc_log<-temp_roc
      roc_log_auc<-temp_roc$auc
      roc_log_cluster<-i
    }
  }
  for(i in levels(mca_2)){
    temp_roc<-roc(Idents(mca_2)==i,GetAssayData(mca_2,'counts')[gene,],quiet = TRUE)
    if((temp_roc$auc>roc_tsnet_auc)&(sum(Idents(mca_2)==i)>5)){
      roc_tsnet<-temp_roc
      roc_tsnet_auc<-temp_roc$auc
      roc_tsnet_cluster<-i
    }
  }
  plot(roc_log,col='red',lwd=3, main=gene,xlab="",ylab="")
  plot(roc_tsnet,col='green',lwd=3,add=TRUE)
  legend(0.77, 0.22, legend=c(paste("Cluster",roc_log_cluster,', AUC:',round(roc_log_auc,2)), paste("Cluster",roc_tsnet_cluster,', AUC:',round(roc_tsnet_auc,2))),
         col=c("red", "green"),lty=1,lwd=3,bty = "n")
  print(roc.test(roc_log,roc_tsnet))
  auck<-as.vector(c(roc_log_auc,roc_tsnet_auc))
  return(auck)
}