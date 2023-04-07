##Load Required Packages
#library('rstudioapi')
library('Seurat')
library('ggplot2')
library('dplyr')
library('pROC')
library('tsne')
library('matlabr')
library('monocle')
library('Nebulosa')
##=========================
## Set default tool paths
tsnetlibpath<-'.'
scentlibpath<-'.'


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
    theme(legend.position = 'none',panel.border = element_rect(linewidth = 2),
          
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(face = 'bold',size=15),
          axis.text.y = element_text(face = 'bold',size=15,colour = 'black'))
  
}