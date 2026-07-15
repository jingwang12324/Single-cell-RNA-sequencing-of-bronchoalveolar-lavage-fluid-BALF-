rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/Myproject/20250421_COPD_BAL_scRNA/")

###scRNA Trans
if(T){
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(gridExtra)

  sce <- read.table('raw_E_matrix.csv',sep = '\t',header = T,row.names = 1,check.names = F)
  counts <- as.matrix(sce)
  colnames(counts) <- colnames(sce)
  rownames(counts) <- rownames(sce)
  write.csv(counts,file = "E_counts.csv")

  sce <- read.table('raw_L_matrix.csv',sep = '\t',header = T,row.names = 1,check.names = F)
  counts <- as.matrix(sce)
  colnames(counts) <- colnames(sce)
  rownames(counts) <- rownames(sce)
  write.csv(counts,file = "L_counts.csv")
  
  sce <- read.table('raw_M_matrix.csv',sep = '\t',header = T,row.names = 1,check.names = F)
  counts <- as.matrix(sce)
  colnames(counts) <- colnames(sce)
  rownames(counts) <- rownames(sce)
  write.csv(counts,file = "M_counts.csv")
  
  sce <- read.table('raw_N_matrix.csv',sep = '\t',header = T,row.names = 1,check.names = F)
  counts <- as.matrix(sce)
  colnames(counts) <- colnames(sce)
  rownames(counts) <- rownames(sce)
  write.csv(counts,file = "N_counts.csv") 
}

###cell singlecell celltype percent
if(T){
  ###Allcell percent
  metadata <- read.table('adata.obs.csv',sep = ',',header = T,check.names = F)
  table(metadata$Group)
  metadata$Celltype
  meta_filt <- metadata[metadata$Group%in%c('COPD','Ctrl'),]
  meta_filt$group <- factor(as.vector(meta_filt$Group),levels=c('COPD','Ctrl'))
  summary <- table(meta_filt[,c('Celltype','Group')])
  df <- data.frame(summary)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  df <- df %>% group_by(Group) %>% mutate(percent = Freq/sum(Freq))
  
  pp <- ggplot(df,aes(x = 3,y = percent,fill= Celltype)) + geom_col(width=1.5,color ='white') + facet_grid(.~Group)
  pp
  colorlist <- c("#dc8e97","#e3d1db","#74a893","#e98741","#5ac6e9","#ebce8e","#7587b1","#c7deef","#e5c06e","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79")
  p1 <- pp + coord_polar(theta = "y") + xlim(c(0.2,3.8)) + 
    scale_fill_manual(values = colorlist) + theme_void() + 
    theme(
      strip.text.x = element_text(size = 14),
      legend.title= element_text(size=15),
      legend.text = element_text(size = 14))
  p1
  ggsave('Allcell_Percent_Circplot.pdf',p1,width = 5,height = 6)
  
  df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
  p <- ggplot(df,aes(x= Group,y=percent,fill=Celltype,stratum=Celltype,alluvium=Celltype)) + 
    scale_fill_manual(values=colorlist) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_classic()
  p5 <- p + geom_col(width = 0.6, color = NA,size = 0.5) + 
    geom_flow(width = 0.6,alpha = 0.22,knot.pos = 0,color = 'white',size = 0.5) + 
    geom_alluvium(width = 0.6,alpha = 1,knot.pos = 0,fill = NA,color = 'white',size = 0.5)
  p5
  ggsave('Allcell_Percent_Barplot.pdf',p5,width = 5,height = 6)
  write.csv(df,'Allcell_Count.csv')
  
  ###Ecell percent
  metadata <- read.table('adata_E.obs.csv',sep = ',',header = T,check.names = F)
  table(metadata$Group)
  metadata$Celltype
  meta_filt <- metadata[metadata$Group%in%c('COPD','Ctrl'),]
  meta_filt$group <- factor(as.vector(meta_filt$Group),levels=c('COPD','Ctrl'))
  summary <- table(meta_filt[,c('Celltype','Group')])
  df <- data.frame(summary)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  df <- df %>% group_by(Group) %>% mutate(percent = Freq/sum(Freq))
  
  pp <- ggplot(df,aes(x = 3,y = percent,fill= Celltype)) + geom_col(width=1.5,color ='white') + facet_grid(.~Group)
  pp
  colorlist <- c("#dc8e97","#e3d1db","#74a893","#5ac6e9","#ebce8e","#e98741","#7587b1","#c7deef","#e5c06e","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79")
  p1 <- pp + coord_polar(theta = "y") + xlim(c(0.2,3.8)) + 
    scale_fill_manual(values = colorlist) + theme_void() + 
    theme(
      strip.text.x = element_text(size = 14),
      legend.title= element_text(size=15),
      legend.text = element_text(size = 14))
  p1
  ggsave('Ecell_Percent_Circplot.pdf',p1,width = 5,height = 6)
  
  df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
  p <- ggplot(df,aes(x= Group,y=percent,fill=Celltype,stratum=Celltype,alluvium=Celltype)) + 
    scale_fill_manual(values=colorlist) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_classic()
  p5 <- p + geom_col(width = 0.6, color = NA,size = 0.5) + 
    geom_flow(width = 0.6,alpha = 0.22,knot.pos = 0,color = 'white',size = 0.5) + 
    geom_alluvium(width = 0.6,alpha = 1,knot.pos = 0,fill = NA,color = 'white',size = 0.5)
  p5
  ggsave('Ecell_Percent_Barplot.pdf',p5,width = 5,height = 6)
  write.csv(df,'Ecell_Group_Count.csv')
  
  ###Lcell percent
  metadata <- read.table('adata_L.obs.csv',sep = ',',header = T,check.names = F)
  table(metadata$Group)
  metadata$Celltype
  meta_filt <- metadata[metadata$Group%in%c('COPD','Ctrl'),]
  meta_filt$group <- factor(as.vector(meta_filt$Group),levels=c('COPD','Ctrl'))
  summary <- table(meta_filt[,c('Celltype','Group')])
  df <- data.frame(summary)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  df <- df %>% group_by(Group) %>% mutate(percent = Freq/sum(Freq))
  
  pp <- ggplot(df,aes(x = 3,y = percent,fill= Celltype)) + geom_col(width=1.5,color ='white') + facet_grid(.~Group)
  pp
  colorlist <- c("#dc8e97","#e3d1db","#74a893","#5ac6e9","#ebce8e","#e98741","#7587b1","#c7deef","#e5c06e","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79")
  p1 <- pp + coord_polar(theta = "y") + xlim(c(0.2,3.8)) + 
    scale_fill_manual(values = colorlist) + theme_void() + 
    theme(
      strip.text.x = element_text(size = 14),
      legend.title= element_text(size=15),
      legend.text = element_text(size = 14))
  p1
  ggsave('Lcell_Percent_Circplot.pdf',p1,width = 5,height = 6)
  
  df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
  p <- ggplot(df,aes(x= Group,y=percent,fill=Celltype,stratum=Celltype,alluvium=Celltype)) + 
    scale_fill_manual(values=colorlist) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_classic()
  p5 <- p + geom_col(width = 0.6, color = NA,size = 0.5) + 
    geom_flow(width = 0.6,alpha = 0.22,knot.pos = 0,color = 'white',size = 0.5) + 
    geom_alluvium(width = 0.6,alpha = 1,knot.pos = 0,fill = NA,color = 'white',size = 0.5)
  p5
  ggsave('Lcell_Percent_Barplot.pdf',p5,width = 5,height = 6)
  write.csv(df,'Lcell_Group_Count.csv') 
  
  ###Mcell percent
  metadata <- read.table('adata_M.obs.csv',sep = ',',header = T,check.names = F)
  table(metadata$Group)
  metadata$Celltype
  meta_filt <- metadata[metadata$Group%in%c('COPD','Ctrl'),]
  meta_filt$group <- factor(as.vector(meta_filt$Group),levels=c('COPD','Ctrl'))
  summary <- table(meta_filt[,c('Celltype','Group')])
  df <- data.frame(summary)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  df <- df %>% group_by(Group) %>% mutate(percent = Freq/sum(Freq))
  
  pp <- ggplot(df,aes(x = 3,y = percent,fill= Celltype)) + geom_col(width=1.5,color ='white') + facet_grid(.~Group)
  colorlist <- c("#dc8e97","#e3d1db","#74a893","#5ac6e9","#ebce8e","#e98741","#7587b1","#c7deef","#e5c06e","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79")
  p1 <- pp + coord_polar(theta = "y") + xlim(c(0.2,3.8)) + 
    scale_fill_manual(values = colorlist) + theme_void() + 
    theme(
      strip.text.x = element_text(size = 14),
      legend.title= element_text(size=15),
      legend.text = element_text(size = 14))
  p1
  ggsave('Mcell_Percent_Circplot.pdf',p1,width = 5,height = 6)
  
  df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
  p <- ggplot(df,aes(x= Group,y=percent,fill=Celltype,stratum=Celltype,alluvium=Celltype)) + 
    scale_fill_manual(values=colorlist) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_classic()
  p5 <- p + geom_col(width = 0.6, color = NA,size = 0.5) + 
    geom_flow(width = 0.6,alpha = 0.22,knot.pos = 0,color = 'white',size = 0.5) + 
    geom_alluvium(width = 0.6,alpha = 1,knot.pos = 0,fill = NA,color = 'white',size = 0.5)
  p5
  ggsave('Mcell_Percent_Barplot.pdf',p5,width = 5,height = 6)
  write.csv(df,'Mcell_Group_Count.csv') 
  
  ###Ncell percent
  metadata <- read.table('adata_N.obs.csv',sep = ',',header = T,check.names = F)
  table(metadata$Group)
  metadata$Celltype
  meta_filt <- metadata[metadata$Group%in%c('COPD','Ctrl'),]
  meta_filt$group <- factor(as.vector(meta_filt$Group),levels=c('COPD','Ctrl'))
  summary <- table(meta_filt[,c('Celltype','Group')])
  df <- data.frame(summary)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  df <- df %>% group_by(Group) %>% mutate(percent = Freq/sum(Freq))
  
  pp <- ggplot(df,aes(x = 3,y = percent,fill= Celltype)) + geom_col(width=1.5,color ='white') + facet_grid(.~Group)
  colorlist <- c("#dc8e97","#e3d1db","#74a893","#5ac6e9","#ebce8e","#e98741","#7587b1","#c7deef","#e5c06e","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79")
  p1 <- pp + coord_polar(theta = "y") + xlim(c(0.2,3.8)) + 
    scale_fill_manual(values = colorlist) + theme_void() + 
    theme(
      strip.text.x = element_text(size = 14),
      legend.title= element_text(size=15),
      legend.text = element_text(size = 14))
  p1
  ggsave('Ncell_Percent_Circplot.pdf',p1,width = 5,height = 6)
  
  df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
  p <- ggplot(df,aes(x= Group,y=percent,fill=Celltype,stratum=Celltype,alluvium=Celltype)) + 
    scale_fill_manual(values=colorlist) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_classic()
  p5 <- p + geom_col(width = 0.6, color = NA,size = 0.5) + 
    geom_flow(width = 0.6,alpha = 0.22,knot.pos = 0,color = 'white',size = 0.5) + 
    geom_alluvium(width = 0.6,alpha = 1,knot.pos = 0,fill = NA,color = 'white',size = 0.5)
  p5
  ggsave('Ncell_Percent_Barplot.pdf',p5,width = 5,height = 6)
  write.csv(df,'Ncell_Group_Count.csv') 
  
  
}

###singlecell Scissor
if(T){
  ###Lsinglecell Scissor
  library(Scissor)
  library(Seurat)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  exprset <- read.table('raw_L_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  metadata <- read.table('adata_L.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  sce <- CreateSeuratObject(exprset,project = 'Lcell',assay = 'RNA',meta.data = metadata)
  sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 1e4)
  sce <- FindVariableFeatures(sce)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- RunUMAP(sce, dims = 1:30)
  sce <- FindNeighbors(sce, dims = 1:10)
  saveRDS(sce,'Lnormalize.sce_v4.rds')

  sce <- readRDS('Lnormalize.sce_v4.rds')
  bulk_exprs <- as.matrix(read.table('GSE30063_expr_ssGSEA.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  pheno <- data.frame(read.table('GSE30063_group.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  
  hvg <- VariableFeatures(sce)
  
  trace("Scissor", edit=TRUE)
  dim(bulk_exprs)
  table(pheno$group)
  scissor_output <- Scissor(bulk_dataset = bulk_exprs, 
                            sc_dataset = sce, 
                            phenotype = pheno$Group1, 
                            tag = c(0, 1), 
                            family = "binomial", 
                            alpha = 0.05, 
                            cutoff = 0.2)
  # 查看Scissor识别的细胞
  Scissor_select <- rep(0, ncol(sce))
  names(Scissor_select) <- colnames(sce)
  Scissor_select[scissor_output$Scissor_pos] <- 1
  Scissor_select[scissor_output$Scissor_neg] <- 2
  sce <- AddMetaData(sce, metadata = Scissor_select, col.name = "scissor")
  
  metadata <- sce@meta.data
  write.csv(metadata,'Lscissor_result.csv')
  saveRDS(scissor_output,'Lscissor_result.rds')

  ###Msinglecell Scissor
  library(Scissor)
  library(Seurat)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  exprset <- read.table('raw_M_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  metadata <- read.table('adata_M.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  sce <- CreateSeuratObject(exprset,project = 'Mcell',assay = 'RNA',meta.data = metadata)
  sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 1e4)
  sce <- FindVariableFeatures(sce)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- RunUMAP(sce, dims = 1:30)
  sce <- FindNeighbors(sce, dims = 1:10)
  saveRDS(sce,'Mnormalize.sce_v4.rds')
  
  sce <- readRDS('Mnormalize.sce_v4.rds')
  bulk_exprs <- as.matrix(read.table('GSE30063_expr_ssGSEA.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  pheno <- data.frame(read.table('GSE30063_group.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  
  #trace("Scissor", edit=TRUE)
  dim(bulk_exprs)
  table(pheno$group)
  scissor_output <- Scissor(bulk_dataset = bulk_exprs, 
                            sc_dataset = sce, 
                            phenotype = pheno$Group1, 
                            tag = c(0, 1), 
                            family = "binomial", 
                            alpha = 0.05, 
                            cutoff = 0.2)
  # 查看Scissor识别的细胞
  Scissor_select <- rep(0, ncol(sce))
  names(Scissor_select) <- colnames(sce)
  Scissor_select[scissor_output$Scissor_pos] <- 1
  Scissor_select[scissor_output$Scissor_neg] <- 2
  sce <- AddMetaData(sce, metadata = Scissor_select, col.name = "scissor")
  
  metadata <- sce@meta.data
  write.csv(metadata,'Mscissor_result.csv')
  saveRDS(scissor_output,'Mscissor_result.rds')
  
  ###Nsinglecell Scissor
  library(Scissor)
  library(Seurat)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  exprset <- read.table('raw_N_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  metadata <- read.table('adata_N.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  sce <- CreateSeuratObject(exprset,project = 'Ncell',assay = 'RNA',meta.data = metadata)
  sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 1e4)
  sce <- FindVariableFeatures(sce)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- RunUMAP(sce, dims = 1:30)
  sce <- FindNeighbors(sce, dims = 1:10)
  saveRDS(sce,'Nnormalize.sce_v4.rds')
  
  sce <- readRDS('Nnormalize.sce_v4.rds')
  bulk_exprs <- as.matrix(read.table('GSE30063_expr_ssGSEA.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  pheno <- data.frame(read.table('GSE30063_group.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  
  #trace("Scissor", edit=TRUE)
  dim(bulk_exprs)
  table(pheno$group)
  scissor_output <- Scissor(bulk_dataset = bulk_exprs, 
                            sc_dataset = sce, 
                            phenotype = pheno$Group1, 
                            tag = c(0, 1), 
                            family = "binomial", 
                            alpha = 0.05, 
                            cutoff = 0.2)
  # 查看Scissor识别的细胞
  Scissor_select <- rep(0, ncol(sce))
  names(Scissor_select) <- colnames(sce)
  Scissor_select[scissor_output$Scissor_pos] <- 1
  Scissor_select[scissor_output$Scissor_neg] <- 2
  sce <- AddMetaData(sce, metadata = Scissor_select, col.name = "scissor")
  
  metadata <- sce@meta.data
  write.csv(metadata,'Nscissor_result.csv')
  saveRDS(scissor_output,'Nscissor_result.rds')
  
  ###Esinglecell Scissor
  library(Scissor)
  library(Seurat)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  exprset <- read.table('raw_E_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  metadata <- read.table('adata_E.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  sce <- CreateSeuratObject(exprset,project = 'Ecell',assay = 'RNA',meta.data = metadata)
  sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 1e4)
  sce <- FindVariableFeatures(sce)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- RunUMAP(sce, dims = 1:30)
  sce <- FindNeighbors(sce, dims = 1:10)
  saveRDS(sce,'Enormalize.sce_v4.rds')
  
  sce <- readRDS('Enormalize.sce_v4.rds')
  bulk_exprs <- as.matrix(read.table('GSE30063_expr_ssGSEA.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  pheno <- data.frame(read.table('GSE30063_group.txt',sep = '\t',header = T,row.names = 1,check.names = F))
  
  #trace("Scissor", edit=TRUE)
  dim(bulk_exprs)
  table(pheno$group)
  scissor_output <- Scissor(bulk_dataset = bulk_exprs, 
                            sc_dataset = sce, 
                            phenotype = pheno$Group1, 
                            tag = c(0, 1), 
                            family = "binomial", 
                            alpha = 0.05, 
                            cutoff = 0.2)
  # 查看Scissor识别的细胞
  Scissor_select <- rep(0, ncol(sce))
  names(Scissor_select) <- colnames(sce)
  Scissor_select[scissor_output$Scissor_pos] <- 1
  Scissor_select[scissor_output$Scissor_neg] <- 2
  sce <- AddMetaData(sce, metadata = Scissor_select, col.name = "scissor")
  
  metadata <- sce@meta.data
  write.csv(metadata,'Escissor_result.csv')
  saveRDS(scissor_output,'Escissor_result.rds')
}

###cellmiloR
if(T){
  ##Ecell miloR
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(zellkonverter)
  library(miloR)
  library(patchwork)
  library(scater)
  sce <- readH5AD('adata_E_result.h5ad')
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_E.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  
  milo.meta <- metadata
  sce
  milo.obj <- Milo(sce)
  milo.obj
  
  milo.obj <- buildGraph(milo.obj, k = 30, d = 30,reduced.dim = "X_pca")
  milo.obj <- makeNhoods(milo.obj, k = 30, d = 30, refined = TRUE, prop = 0.05,reduced_dims = "X_pca")
  plotNhoodSizeHist(milo.obj)
  #trace(miloR:::calcNhoodDistance, edit = T)
  #milo.obj <- calcNhoodDistance(milo.obj, d=30,reduced.dim = "X_pca")
  milo.obj <- countCells(milo.obj, samples="batch", meta.data=milo.meta)
  
  milo.design <- as.data.frame(xtabs(~ Group + batch, data=milo.meta))
  milo.design <- milo.design[milo.design$Freq > 0, ]
  rownames(milo.design) <- milo.design$batch
  milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
  
  milo.res <- testNhoods(milo.obj, design=~Group, design.df=milo.design,reduced.dim = "X_pca")
  head(milo.res)
  write.table(milo.res,'Emilo.res.xls',sep = '\t')
  saveRDS(milo.obj,'Emilo.obj.rds')
  
  milo.res %>%arrange(SpatialFDR) %>%head() 
  ggplot(milo.res, aes(PValue)) + geom_histogram(bins=50)
  ggplot(milo.res, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
  
  #milo.obj <- readRDS('Emilo.obj.rds')
  #milo.res <- read.table('Emilo.res.xls',sep = '\t')
  milo.obj <- buildNhoodGraph(milo.obj)
  ## Plot single-cell UMAP
  milo.obj
  umap_pl <- plotReducedDim(milo.obj, dimred = "X_umap", 
                            colour_by = "Group", text_by = "Celltype", 
                            text_size = 3, point_size=0.5) + guides(fill="none")
  ##umap_pl
  ## Plot neighbourhood graph
  nh_graph_pl <- plotNhoodGraphDA(milo.obj, milo.res,layout="X_umap",alpha=0.1)  #alpha默认0.1
  p1 <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
  p1
  ggsave('Emilo_Umap_Result.pdf',p1,width = 10,height = 6)
  da_results <- annotateNhoods(milo.obj, milo.res, coldata_col ="Celltype")
  
  ggplot(da_results, aes(Celltype_fraction)) + geom_histogram(bins=50)
  da_results$Celltype <- ifelse(da_results$Celltype_fraction <0.7,"Mixed", da_results$Celltype)
  p1 <- plotDAbeeswarm(da_results, group.by = "Celltype")
  ggsave('EDAbeeswarm.pdf',p1,width = 6,height = 5)

  ##Lcell miloR
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(zellkonverter)
  library(miloR)
  library(patchwork)
  library(scater)
  sce <- readH5AD('adata_L_result.h5ad')
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_L.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  
  milo.meta <- metadata
  sce
  milo.obj <- Milo(sce)
  milo.obj
  
  milo.obj <- buildGraph(milo.obj, k = 30, d = 30,reduced.dim = "X_pca")
  milo.obj <- makeNhoods(milo.obj, k = 30, d = 30, refined = TRUE, prop = 0.05,reduced_dims = "X_pca")
  plotNhoodSizeHist(milo.obj)
  #trace(miloR:::calcNhoodDistance, edit = T)
  #milo.obj <- calcNhoodDistance(milo.obj, d=30,reduced.dim = "X_pca")
  milo.obj <- countCells(milo.obj, samples="batch", meta.data=milo.meta)
  
  milo.design <- as.data.frame(xtabs(~ Group + batch, data=milo.meta))
  milo.design <- milo.design[milo.design$Freq > 0, ]
  rownames(milo.design) <- milo.design$batch
  milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
  
  milo.res <- testNhoods(milo.obj, design=~Group, design.df=milo.design,reduced.dim = "X_pca")
  head(milo.res)
  write.table(milo.res,'Lmilo.res.xls',sep = '\t')
  saveRDS(milo.obj,'Lmilo.obj.rds')
  
  milo.res %>%arrange(SpatialFDR) %>%head() 
  ggplot(milo.res, aes(PValue)) + geom_histogram(bins=50)
  ggplot(milo.res, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
  
  #milo.obj <- readRDS('Emilo.obj.rds')
  #milo.res <- read.table('Emilo.res.xls',sep = '\t')
  milo.obj <- buildNhoodGraph(milo.obj)
  ## Plot single-cell UMAP
  milo.obj
  umap_pl <- plotReducedDim(milo.obj, dimred = "X_umap", 
                            colour_by = "Group", text_by = "Celltype", 
                            text_size = 3, point_size=0.5) + guides(fill="none")
  ##umap_pl
  ## Plot neighbourhood graph
  nh_graph_pl <- plotNhoodGraphDA(milo.obj, milo.res,layout="X_umap",alpha=0.1)  #alpha默认0.1
  p1 <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
  p1
  ggsave('Lmilo_Umap_Result.pdf',p1,width = 10,height = 6)
  da_results <- annotateNhoods(milo.obj, milo.res, coldata_col ="Celltype")
  
  ggplot(da_results, aes(Celltype_fraction)) + geom_histogram(bins=50)
  da_results$Celltype <- ifelse(da_results$Celltype_fraction <0.7,"Mixed", da_results$Celltype)
  p1 <- plotDAbeeswarm(da_results, group.by = "Celltype")
  ggsave('LDAbeeswarm.pdf',p1,width = 6,height = 5)
  
  ##Mcell miloR
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(zellkonverter)
  library(miloR)
  library(patchwork)
  library(scater)
  sce <- readH5AD('adata_M_result.h5ad')
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_M.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  
  milo.meta <- metadata
  sce
  milo.obj <- Milo(sce)
  milo.obj
  
  milo.obj <- buildGraph(milo.obj, k = 30, d = 30,reduced.dim = "X_pca")
  milo.obj <- makeNhoods(milo.obj, k = 30, d = 30, refined = TRUE, prop = 0.05,reduced_dims = "X_pca")
  plotNhoodSizeHist(milo.obj)
  #trace(miloR:::calcNhoodDistance, edit = T)
  #milo.obj <- calcNhoodDistance(milo.obj, d=30,reduced.dim = "X_pca")
  milo.obj <- countCells(milo.obj, samples="batch", meta.data=milo.meta)
  
  milo.design <- as.data.frame(xtabs(~ Group + batch, data=milo.meta))
  milo.design <- milo.design[milo.design$Freq > 0, ]
  rownames(milo.design) <- milo.design$batch
  milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
  
  milo.res <- testNhoods(milo.obj, design=~Group, design.df=milo.design,reduced.dim = "X_pca")
  head(milo.res)
  write.table(milo.res,'Mmilo.res.xls',sep = '\t')
  saveRDS(milo.obj,'Mmilo.obj.rds')
  
  milo.res %>%arrange(SpatialFDR) %>%head() 
  ggplot(milo.res, aes(PValue)) + geom_histogram(bins=50)
  ggplot(milo.res, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
  
  #milo.obj <- readRDS('Emilo.obj.rds')
  #milo.res <- read.table('Emilo.res.xls',sep = '\t')
  milo.obj <- buildNhoodGraph(milo.obj)
  ## Plot single-cell UMAP
  milo.obj
  umap_pl <- plotReducedDim(milo.obj, dimred = "X_umap", 
                            colour_by = "Group", text_by = "Celltype", 
                            text_size = 3, point_size=0.5) + guides(fill="none")
  ##umap_pl
  ## Plot neighbourhood graph
  nh_graph_pl <- plotNhoodGraphDA(milo.obj, milo.res,layout="X_umap",alpha=0.1)  #alpha默认0.1
  p1 <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
  p1
  ggsave('Mmilo_Umap_Result.pdf',p1,width = 10,height = 6)
  da_results <- annotateNhoods(milo.obj, milo.res, coldata_col ="Celltype")
  
  ggplot(da_results, aes(Celltype_fraction)) + geom_histogram(bins=50)
  da_results$Celltype <- ifelse(da_results$Celltype_fraction <0.7,"Mixed", da_results$Celltype)
  p1 <- plotDAbeeswarm(da_results, group.by = "Celltype")
  ggsave('MDAbeeswarm.pdf',p1,width = 6,height = 5)
  
  ##Ncell miloR
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(zellkonverter)
  library(miloR)
  library(patchwork)
  library(scater)
  sce <- readH5AD('adata_N_result.h5ad')
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_N.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  
  milo.meta <- metadata
  sce
  milo.obj <- Milo(sce)
  milo.obj
  
  milo.obj <- buildGraph(milo.obj, k = 30, d = 30,reduced.dim = "X_pca")
  milo.obj <- makeNhoods(milo.obj, k = 30, d = 30, refined = TRUE, prop = 0.05,reduced_dims = "X_pca")
  plotNhoodSizeHist(milo.obj)
  #trace(miloR:::calcNhoodDistance, edit = T)
  #milo.obj <- calcNhoodDistance(milo.obj, d=30,reduced.dim = "X_pca")
  milo.obj <- countCells(milo.obj, samples="batch", meta.data=milo.meta)
  
  milo.design <- as.data.frame(xtabs(~ Group + batch, data=milo.meta))
  milo.design <- milo.design[milo.design$Freq > 0, ]
  rownames(milo.design) <- milo.design$batch
  milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
  
  milo.res <- testNhoods(milo.obj, design=~Group, design.df=milo.design,reduced.dim = "X_pca")
  head(milo.res)
  write.table(milo.res,'Nmilo.res.xls',sep = '\t')
  saveRDS(milo.obj,'Nmilo.obj.rds')
  
  milo.res %>%arrange(SpatialFDR) %>%head() 
  ggplot(milo.res, aes(PValue)) + geom_histogram(bins=50)
  ggplot(milo.res, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
  
  #milo.obj <- readRDS('Emilo.obj.rds')
  #milo.res <- read.table('Emilo.res.xls',sep = '\t')
  milo.obj <- buildNhoodGraph(milo.obj)
  ## Plot single-cell UMAP
  milo.obj
  umap_pl <- plotReducedDim(milo.obj, dimred = "X_umap", 
                            colour_by = "Group", text_by = "Celltype", 
                            text_size = 3, point_size=0.5) + guides(fill="none")
  ##umap_pl
  ## Plot neighbourhood graph
  nh_graph_pl <- plotNhoodGraphDA(milo.obj, milo.res,layout="X_umap",alpha=0.1)  #alpha默认0.1
  p1 <- umap_pl + nh_graph_pl + plot_layout(guides="collect")
  p1
  ggsave('Nmilo_Umap_Result.pdf',p1,width = 10,height = 6)
  da_results <- annotateNhoods(milo.obj, milo.res, coldata_col ="Celltype")
  
  ggplot(da_results, aes(Celltype_fraction)) + geom_histogram(bins=50)
  da_results$Celltype <- ifelse(da_results$Celltype_fraction <0.7,"Mixed", da_results$Celltype)
  p1 <- plotDAbeeswarm(da_results, group.by = "Celltype")
  p1
  ggsave('NDAbeeswarm.pdf',p1,width = 6,height = 5)
  }

###SCENIC
if(T){
  ###Epithelial SCENIC
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(gridExtra)
  library(zellkonverter)
  exprset <- read.table('raw_E_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_E.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  E_sce <- CreateSeuratObject(exprset,project = 'E',assay = 'RNA',meta.data = metadata)
  E_sce@meta.data$Celltype <- metadata$Celltype2
  loom <- open_loom("SCENIC/E_aucell.loom")
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  #将调节因子矩阵转换为基因列表
  regulons <- regulonsToGeneLists(regulons_incidMat)
  #获取调节因子的AUC值
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  #获取调节因子的阈值
  regulonAucThresholds <- get_regulon_thresholds(loom)
  sub_regulonAUC <- regulonAUC[,match(colnames(E_sce),colnames(regulonAUC))]
  identical(colnames(sub_regulonAUC), colnames(E_sce))
  #每个细胞的aucell值添加到seurat对象的meta.data中
  E_sce@meta.data = cbind(E_sce@meta.data ,t(sub_regulonAUC@assays@data$AUC))
  
  rss<-calcRSS(AUC=getAUC(sub_regulonAUC),cellAnnotation=E_sce@meta.data$Celltype)
  rssPlot <- plotRSS(rss)
  rssPlot$plot
  rss_data = rss
  important_tfs <- lapply(colnames(rss_data), function(celltype){
    head(rownames(rss_data)[order(rss_data[, celltype], decreasing = TRUE)], 5)
  })
  important_tfs <- unique(unlist(important_tfs))
  
  rssPlot <- plotRSS(rss[important_tfs, ])
  # 进一步美化图形
  optimized_plot <- rssPlot$plot + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Transcription Factors by RSS") +
    scale_fill_viridis_c(option = "plasma")  # 使用viridis配色方案
  
  # 显示图形
  optimized_plot
  ggsave('E_SCENIC_Transcription Factors by RSS Dotplot.pdf',optimized_plot,width = 5,height = 6)
  table(E_sce@meta.data$Celltype)
  #分细胞类型展示regulon特异性分数
  p1 <- plotRSS_oneSet(rss,'Basal cell')
  p2 <- plotRSS_oneSet(rss,'Ciliated cell')
  p3 <- plotRSS_oneSet(rss,'Club cell')
  p4 <- plotRSS_oneSet(rss,'Epithelial cell')
  p5 <- plotRSS_oneSet(rss,'Fibroblast')
  p6 <- plotRSS_oneSet(rss,'Mesenchymal cell')
  p1 <- grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 4, nrow = 2)
  p1
  #ggsave('E_SCENIC_RSS_Transcription Factors.pdf',p1,width = 10,height = 11)
  regulonAUC <- getAUC(sub_regulonAUC)
  n_top_regulons <- 50
  regulonAUC <- regulonAUC[names(sort(apply(regulonAUC, 1, sd), decreasing = TRUE))[1:n_top_regulons], ]
  
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),  # Z-score 范围
    colors = c("white", "gray90", "black") 
  )
  # 绘制热图
  acell_scaled_data <- t(scale(t(regulonAUC)))
  pdf('E_SCENIC_Allcell_Heatmap.pdf',width = 5,height = 5)
  pheatmap(acell_scaled_data, 
           scale = "row", 
           clustering_method = "ward.D2",
           color = col_fun,
           show_colnames = FALSE,
           fontsize_row = 6)
  dev.off()
  
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_types <- E_sce$Celltype # 假设seurat对象中有celltype信息
  
  regulon_activity_by_celltype <- sapply(split(colnames(regulonAUC), cell_types), 
                                         function(cells) rowMeans(regulonAUC[,cells]))
  
  # 选择最具细胞类型特异性的regulon
  top_regulon_per_celltype <- apply(regulon_activity_by_celltype, 2, 
                                    function(x) names(x)[which.max(x)])
  cluster_to_celltype <- c(
    "Basal cell" = "Basal cell",
    "Ciliated cell" = "Ciliated cell",
    "Club cell" = "Club cell",
    "Epithelial cell" = "Epithelial cell",
    "Fibroblast" = "Fibroblast",
    "Mesenchymal cell" = "Mesenchymal cell"
  )
  ###"#dc8e97","#e3d1db","#74a893","#e98741","#5ac6e9","#ebce8e","#e5c06e","#7587b1","#c7deef","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79"
  celltype_colors<- c(
    "Basal cell" = "#dc8e97",
    "Ciliated cell" = "#e3d1db",
    "Club cell" = "#e5c06e",
    "Epithelial cell" = "#74a893",
    "Fibroblast" = "#e98741",
    "Mesenchymal cell" = "#5ac6e9"
  )
  # 假设 heatmap_columns 是热图的列名（聚类名称）
  heatmap_columns <- colnames(regulon_activity_by_celltype)
  
  # 创建列注释
  column_ha <- HeatmapAnnotation(
    CellType = cluster_to_celltype[heatmap_columns],
    col = list(CellType = celltype_colors),  # 颜色映射
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  scaled_data <- t(scale(t(regulon_activity_by_celltype)))
  
  # 绘制热图
  pdf('E_SCENIC_Celltype_Heatmap.pdf',width = 6,height = 7)
  Heatmap(scaled_data,
          name = "Regulon activity",
          top_annotation = column_ha,
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_Group <- E_sce$Group # 假设seurat对象中有celltype信息
  regulon_activity_by_Group <- sapply(split(colnames(regulonAUC), cell_Group), 
                                      function(cells) rowMeans(regulonAUC[,cells]))
  # 选择最具细胞类型特异性的regulon
  top_regulon_Group  <- apply(regulon_activity_by_Group, 2, 
                              function(x) names(x)[which.max(x)])
  #scaled_data_Group <- t(scale(t(regulon_activity_by_Group)))
  write.table(regulon_activity_by_Group,'E_SCENIC_Group.xls',sep = '\t') 
  pdf('E_SCENIC_Group_Heatmap.pdf',width = 6,height = 6)
  Heatmap(regulon_activity_by_Group,
          name = "Regulon activity",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  write.table(E_sce@meta.data,'E_SCENIC_sce.xls',sep = '\t') 
  
  grn <- read.table('T_grn_output.tsv',sep = '\t',check.names = F,header = T) 
  pparg <- grn[which(grn$TF == 'PPARG'),]
  elf1 <- grn[which(grn$TF == 'ELF1'),]
  write.table(pparg,'TH1_PPARG_TF_Target.xls',sep = '\t')
  write.table(elf1,'TH1_ELF1_TF_Target.xls',sep = '\t')
  
  library(ggplot2)
  library(ggridges) # 用于绘制山峦图
  library(dplyr)    # 用于数据处理
  library(tidyr)    # 用于数据整理
  data <- read.table('E_SCENIC_sce.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  data <- data[,c(-1:-3,-5:-13)]
  data_long <- data %>%
    pivot_longer(cols = -Group, names_to = "TF", values_to = "AUC")
  
  ###RFX3(+) BCL6(+)
  selected_TF <- "RFX3(+)"
  pdf('E_Ridges_Activity of RFX3(+) Across Group.pdf',width = 5,height = 4)
   data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()

  selected_TF <- "BCL6(+)"
  pdf('E_Ridges_Activity of BCL6(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
  
  ###Monocytes SCENIC
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(gridExtra)
  library(zellkonverter)
  exprset <- read.table('raw_M_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_M.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  M_sce <- CreateSeuratObject(exprset,project = 'M',assay = 'RNA',meta.data = metadata)

  loom <- open_loom("SCENIC/M_aucell.loom")
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  #将调节因子矩阵转换为基因列表
  regulons <- regulonsToGeneLists(regulons_incidMat)
  #获取调节因子的AUC值
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  #获取调节因子的阈值
  regulonAucThresholds <- get_regulon_thresholds(loom)
  sub_regulonAUC <- regulonAUC[,match(colnames(M_sce),colnames(regulonAUC))]
  identical(colnames(sub_regulonAUC), colnames(M_sce))
  #每个细胞的aucell值添加到seurat对象的meta.data中
  M_sce@meta.data = cbind(M_sce@meta.data ,t(sub_regulonAUC@assays@data$AUC))
  
  rss<-calcRSS(AUC=getAUC(sub_regulonAUC),cellAnnotation=M_sce@meta.data$Celltype)
  rssPlot <- plotRSS(rss)
  rssPlot$plot
  rss_data = rss
  important_tfs <- lapply(colnames(rss_data), function(celltype){
    head(rownames(rss_data)[order(rss_data[, celltype], decreasing = TRUE)], 5)
  })
  important_tfs <- unique(unlist(important_tfs))
  
  rssPlot <- plotRSS(rss[important_tfs, ])
  # 进一步美化图形
  optimized_plot <- rssPlot$plot + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Transcription Factors by RSS") +
    scale_fill_viridis_c(option = "plasma")  # 使用viridis配色方案
  optimized_plot
  ggsave('M_SCENIC_Transcription Factors by RSS Dotplot.pdf',optimized_plot,width = 5,height = 6)
  
  table(M_sce@meta.data$Celltype)
  #分细胞类型展示regulon特异性分数
  p1 <- plotRSS_oneSet(rss,'AMs')
  p2 <- plotRSS_oneSet(rss,'cDC1')
  p3 <- plotRSS_oneSet(rss,'cDC2')
  p4 <- plotRSS_oneSet(rss,'Macro-C1QC')
  p5 <- plotRSS_oneSet(rss,'Macro-CD163')
  p6 <- plotRSS_oneSet(rss,'Macro-FABP4')
  p7 <- plotRSS_oneSet(rss,'Macro-FCN1')
  p8 <- plotRSS_oneSet(rss,'Macro-G0S2')
  p9 <- plotRSS_oneSet(rss,'Macro-IL1B')
  p10 <- plotRSS_oneSet(rss,'Macro-MARCO')
  p11 <- plotRSS_oneSet(rss,'Macro-SPP1')
  p12 <- plotRSS_oneSet(rss,'Macro-TNFSF10')
  p13 <- plotRSS_oneSet(rss,'Mast')
  p14 <- plotRSS_oneSet(rss,'pDC')
  p1 <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4, nrow = 2)
  p1
  p2 <- grid.arrange(p9,p10,p11,p12,p13,p14, ncol = 4, nrow = 2)
  p2
  #ggsave('M_SCENIC_RSS_Transcription Factors.pdf',p1,width = 10,height = 11)
  regulonAUC <- getAUC(sub_regulonAUC)
  n_top_regulons <- 50
  regulonAUC <- regulonAUC[names(sort(apply(regulonAUC, 1, sd), decreasing = TRUE))[1:n_top_regulons], ]
  
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),  # Z-score 范围
    colors = c("white", "gray90", "black") 
  )
  # 绘制热图
  acell_scaled_data <- t(scale(t(regulonAUC)))
  pdf('M_SCENIC_Allcell_Heatmap.pdf',width = 5,height = 5)
  pheatmap(acell_scaled_data, 
           scale = "row", 
           clustering_method = "ward.D2",
           color = col_fun,
           show_colnames = FALSE,
           fontsize_row = 6)
  dev.off()
  
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_types <- M_sce$Celltype # 假设seurat对象中有celltype信息
  
  regulon_activity_by_celltype <- sapply(split(colnames(regulonAUC), cell_types), 
                                         function(cells) rowMeans(regulonAUC[,cells]))
  
  # 选择最具细胞类型特异性的regulon
  top_regulon_per_celltype <- apply(regulon_activity_by_celltype, 2, 
                                    function(x) names(x)[which.max(x)])
  cluster_to_celltype <- c(
    "AMs" = "AMs",
    "cDC1" = "cDC1",
    "cDC2" = "cDC2",
    "Macro-C1QC" = "Macro-C1QC",
    "Macro-CD163" = "Macro-CD163",
    "Macro-FABP4" = "Macro-FABP4",
    "Macro-FCN1" = "Macro-FCN1",
    "Macro-G0S2" = "Macro-G0S2",
    "Macro-IL1B" = "Macro-IL1B",
    "Macro-SPP1" = "Macro-SPP1",
    "Macro-TNFSF10" = "Macro-TNFSF10",
    "Macro-FABP4" = "Macro-FABP4",
    "Mast" = "Mast",
    "pDC" = "pDC"
  )
  ###"#dc8e97","#e3d1db","#74a893","#e98741","#5ac6e9","#ebce8e","#e5c06e","#7587b1","#c7deef","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79"
  celltype_colors<- c(
    "AMs" = "#dc8e97",
    "cDC1" = "#e3d1db",
    "cDC2" = "#74a893",
    "Macro-C1QC" = "#e98741",
    "Macro-CD163" = "#5ac6e9",
    "Macro-FABP4" = "#ebce8e",
    "Macro-FCN1" = "#e5c06e",
    "Macro-G0S2" = "#7587b1",
    "Macro-IL1B" = "#c7deef",
    "Macro-SPP1" = "#e97371",
    "Macro-TNFSF10" = "#e1a4c6",
    "Macro-FABP4" = "#916ba6",
    "Mast" = "#cb8f82",
    "pDC" = "#7db3af"
  )
  # 假设 heatmap_columns 是热图的列名（聚类名称）
  heatmap_columns <- colnames(regulon_activity_by_celltype)
  
  # 创建列注释
  column_ha <- HeatmapAnnotation(
    CellType = cluster_to_celltype[heatmap_columns],
    col = list(CellType = celltype_colors),  # 颜色映射
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  scaled_data <- t(scale(t(regulon_activity_by_celltype)))
  
  # 绘制热图
  pdf('M_SCENIC_Celltype_Heatmap.pdf',width = 7,height = 7)
  Heatmap(scaled_data,
          name = "Regulon activity",
          top_annotation = column_ha,
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_Group <- M_sce$Group # 假设seurat对象中有celltype信息
  regulon_activity_by_Group <- sapply(split(colnames(regulonAUC), cell_Group), 
                                      function(cells) rowMeans(regulonAUC[,cells]))
  # 选择最具细胞类型特异性的regulon
  top_regulon_Group  <- apply(regulon_activity_by_Group, 2, 
                              function(x) names(x)[which.max(x)])
  #scaled_data_Group <- t(scale(t(regulon_activity_by_Group)))
  write.table(regulon_activity_by_Group,'M_SCENIC_Group.xls',sep = '\t') 
  pdf('M_SCENIC_Group_Heatmap.pdf',width = 6,height = 6)
  Heatmap(regulon_activity_by_Group,
          name = "Regulon activity",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  write.table(M_sce@meta.data,'M_SCENIC_sce.xls',sep = '\t') 
  
  grn <- read.table('T_grn_output.tsv',sep = '\t',check.names = F,header = T) 
  pparg <- grn[which(grn$TF == 'PPARG'),]
  elf1 <- grn[which(grn$TF == 'ELF1'),]
  write.table(pparg,'TH1_PPARG_TF_Target.xls',sep = '\t')
  write.table(elf1,'TH1_ELF1_TF_Target.xls',sep = '\t')
  
  library(ggplot2)
  library(ggridges) # 用于绘制山峦图
  library(dplyr)    # 用于数据处理
  library(tidyr)    # 用于数据整理
  data <- read.table('M_SCENIC_sce.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  data <- data[,c(-1:-3,-5:-13)]
  data_long <- data %>%
    pivot_longer(cols = -Group, names_to = "TF", values_to = "AUC")
  
  ###CEBPB(+) CHURC1(+)
  selected_TF <- "CEBPB(+)"
  pdf('M_Ridges_Activity of CEBPB(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
  
  selected_TF <- "CHURC1(+)"
  pdf('M_Ridges_Activity of CHURC1(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
  
  ###Lymphocyte SCENIC
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(gridExtra)
  library(zellkonverter)
  exprset <- read.table('raw_L_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_L.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  L_sce <- CreateSeuratObject(exprset,project = 'L',assay = 'RNA',meta.data = metadata)
  
  loom <- open_loom("SCENIC/L_aucell.loom")
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  #将调节因子矩阵转换为基因列表
  regulons <- regulonsToGeneLists(regulons_incidMat)
  #获取调节因子的AUC值
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  #获取调节因子的阈值
  regulonAucThresholds <- get_regulon_thresholds(loom)
  sub_regulonAUC <- regulonAUC[,match(colnames(L_sce),colnames(regulonAUC))]
  identical(colnames(sub_regulonAUC), colnames(L_sce))
  #每个细胞的aucell值添加到seurat对象的meta.data中
  L_sce@meta.data = cbind(L_sce@meta.data ,t(sub_regulonAUC@assays@data$AUC))
  
  rss<-calcRSS(AUC=getAUC(sub_regulonAUC),cellAnnotation=L_sce@meta.data$Celltype)
  rssPlot <- plotRSS(rss)
  rssPlot$plot
  rss_data = rss
  important_tfs <- lapply(colnames(rss_data), function(celltype){
    head(rownames(rss_data)[order(rss_data[, celltype], decreasing = TRUE)], 5)
  })
  important_tfs <- unique(unlist(important_tfs))
  
  rssPlot <- plotRSS(rss[important_tfs, ])
  # 进一步美化图形
  optimized_plot <- rssPlot$plot + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Transcription Factors by RSS") +
    scale_fill_viridis_c(option = "plasma")  # 使用viridis配色方案
  optimized_plot
  ggsave('L_SCENIC_Transcription Factors by RSS Dotplot.pdf',optimized_plot,width = 5,height = 6)
  
  table(L_sce@meta.data$Celltype)
  #分细胞类型展示regulon特异性分数
  p1 <- plotRSS_oneSet(rss,'B cell')
  p2 <- plotRSS_oneSet(rss,'CD4Tn')
  p3 <- plotRSS_oneSet(rss,'CD8Tcm')
  p4 <- plotRSS_oneSet(rss,'CD8Teff')
  p5 <- plotRSS_oneSet(rss,'CD8Tem')
  p6 <- plotRSS_oneSet(rss,'CD8Tex')
  p7 <- plotRSS_oneSet(rss,'CD8Trm')
  p8 <- plotRSS_oneSet(rss,'MAIT')
  p9 <- plotRSS_oneSet(rss,'NK')
  p10 <- plotRSS_oneSet(rss,'NKT')
  
  p1 <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4, nrow = 2)
  p1
  p2 <- grid.arrange(p9,p10, ncol = 4, nrow = 2)
  p2
  #ggsave('L_SCENIC_RSS_Transcription Factors.pdf',p1,width = 10,height = 11)
  regulonAUC <- getAUC(sub_regulonAUC)
  n_top_regulons <- 50
  regulonAUC <- regulonAUC[names(sort(apply(regulonAUC, 1, sd), decreasing = TRUE))[1:n_top_regulons], ]
  
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),  # Z-score 范围
    colors = c("white", "gray90", "black") 
  )
  # 绘制热图
  acell_scaled_data <- t(scale(t(regulonAUC)))
  pdf('L_SCENIC_Allcell_Heatmap.pdf',width = 5,height = 5)
  pheatmap(acell_scaled_data, 
           scale = "row", 
           clustering_method = "ward.D2",
           color = col_fun,
           show_colnames = FALSE,
           fontsize_row = 6)
  dev.off()
  
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_types <- L_sce$Celltype # 假设seurat对象中有celltype信息
  
  regulon_activity_by_celltype <- sapply(split(colnames(regulonAUC), cell_types), 
                                         function(cells) rowMeans(regulonAUC[,cells]))
  
  # 选择最具细胞类型特异性的regulon
  top_regulon_per_celltype <- apply(regulon_activity_by_celltype, 2, 
                                    function(x) names(x)[which.max(x)])
  cluster_to_celltype <- c(
    "B cell" = "B cell",
    "CD4Tn" = "CD4Tn",
    "CD8Tcm" = "CD8Tcm",
    "CD8Teff" = "CD8Teff",
    "CD8Tem" = "CD8Tem",
    "CD8Tex" = "CD8Tex",
    "CD8Trm" = "CD8Trm",
    "MAIT" = "MAIT",
    "NK" = "NK",
    "NKT" = "NKT"
  )
  ###"#dc8e97","#e3d1db","#74a893","#e98741","#5ac6e9","#ebce8e","#e5c06e","#7587b1","#c7deef","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79"
  celltype_colors<- c(
    "B cell" = "#dc8e97",
    "CD4Tn" = "#e3d1db",
    "CD8Tcm" = "#74a893",
    "CD8Teff" = "#e98741",
    "CD8Tem" = "#5ac6e9",
    "CD8Tex" = "#ebce8e",
    "CD8Trm" = "#e5c06e",
    "MAIT" = "#7587b1",
    "NK" = "#c7deef",
    "NKT" = "#e97371"
  )
  # 假设 heatmap_columns 是热图的列名（聚类名称）
  heatmap_columns <- colnames(regulon_activity_by_celltype)
  
  # 创建列注释
  column_ha <- HeatmapAnnotation(
    CellType = cluster_to_celltype[heatmap_columns],
    col = list(CellType = celltype_colors),  # 颜色映射
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  scaled_data <- t(scale(t(regulon_activity_by_celltype)))
  
  # 绘制热图
  pdf('L_SCENIC_Celltype_Heatmap.pdf',width = 7,height = 7)
  Heatmap(scaled_data,
          name = "Regulon activity",
          top_annotation = column_ha,
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_Group <- L_sce$Group # 假设seurat对象中有celltype信息
  regulon_activity_by_Group <- sapply(split(colnames(regulonAUC), cell_Group), 
                                      function(cells) rowMeans(regulonAUC[,cells]))
  # 选择最具细胞类型特异性的regulon
  top_regulon_Group  <- apply(regulon_activity_by_Group, 2, 
                              function(x) names(x)[which.max(x)])
  #scaled_data_Group <- t(scale(t(regulon_activity_by_Group)))
  write.table(regulon_activity_by_Group,'L_SCENIC_Group.xls',sep = '\t') 
  pdf('L_SCENIC_Group_Heatmap.pdf',width = 6,height = 6)
  Heatmap(regulon_activity_by_Group,
          name = "Regulon activity",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  write.table(L_sce@meta.data,'L_SCENIC_sce.xls',sep = '\t') 
  
  grn <- read.table('T_grn_output.tsv',sep = '\t',check.names = F,header = T) 
  pparg <- grn[which(grn$TF == 'PPARG'),]
  elf1 <- grn[which(grn$TF == 'ELF1'),]
  write.table(pparg,'TH1_PPARG_TF_Target.xls',sep = '\t')
  write.table(elf1,'TH1_ELF1_TF_Target.xls',sep = '\t')
  
  library(ggplot2)
  library(ggridges) # 用于绘制山峦图
  library(dplyr)    # 用于数据处理
  library(tidyr)    # 用于数据整理
  data <- read.table('L_SCENIC_sce.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  data <- data[,c(-1:-3,-5:-12)]
  data_long <- data %>%
    pivot_longer(cols = -Group, names_to = "TF", values_to = "AUC")
  
  ###TAF6(+) BCL6(+)
  selected_TF <- "TAF6(+)"
  pdf('L_Ridges_Activity of TAF6(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
  
  selected_TF <- "BCL6(+)"
  pdf('L_Ridges_Activity of BCL6(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
  
  ###Granulocytes SCENIC
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(biomaRt)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(pheatmap)
  library(circlize)
  library(gridExtra)
  library(zellkonverter)
  exprset <- read.table('raw_N_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  #example_h5ad_milo <- Milo(sce)
  metadata <- read.table('adata_N.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  G_sce <- CreateSeuratObject(exprset,project = 'G',assay = 'RNA',meta.data = metadata)
  
  loom <- open_loom("SCENIC/N_aucell.loom")
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  #将调节因子矩阵转换为基因列表
  regulons <- regulonsToGeneLists(regulons_incidMat)
  #获取调节因子的AUC值
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  #获取调节因子的阈值
  regulonAucThresholds <- get_regulon_thresholds(loom)
  sub_regulonAUC <- regulonAUC[,match(colnames(G_sce),colnames(regulonAUC))]
  identical(colnames(sub_regulonAUC), colnames(G_sce))
  #每个细胞的aucell值添加到seurat对象的meta.data中
  G_sce@meta.data = cbind(G_sce@meta.data ,t(sub_regulonAUC@assays@data$AUC))
  
  rss<-calcRSS(AUC=getAUC(sub_regulonAUC),cellAnnotation=G_sce@meta.data$Celltype)
  rssPlot <- plotRSS(rss)
  rssPlot$plot
  rss_data = rss
  important_tfs <- lapply(colnames(rss_data), function(celltype){
    head(rownames(rss_data)[order(rss_data[, celltype], decreasing = TRUE)], 5)
  })
  important_tfs <- unique(unlist(important_tfs))
  
  rssPlot <- plotRSS(rss[important_tfs, ])
  # 进一步美化图形
  optimized_plot <- rssPlot$plot + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle("Transcription Factors by RSS") +
    scale_fill_viridis_c(option = "plasma")  # 使用viridis配色方案
  optimized_plot
  ggsave('G_SCENIC_Transcription Factors by RSS Dotplot.pdf',optimized_plot,width = 5,height = 6)
  
  table(G_sce@meta.data$Celltype)
  write.table(rss,'Grss.xls',sep = '\t')
  rss <- data.matrix(read.table('Grss.xls',sep = '\t',header = T,row.names = 1,check.names = F))
  #分细胞类型展示regulon特异性分数
  p1 <- plotRSS_oneSet(rss,'Basophil')
  p2 <- plotRSS_oneSet(rss,'Eosinophil')
  p3 <- plotRSS_oneSet(rss,'Neutrophil')
  
  p1 <- grid.arrange(p1,p2,p3, ncol = 4, nrow = 2)
  p1

  #ggsave('G_SCENIC_RSS_Transcription Factors.pdf',p1,width = 10,height = 11)
  regulonAUC <- getAUC(sub_regulonAUC)
  n_top_regulons <- 50
  regulonAUC <- regulonAUC[names(sort(apply(regulonAUC, 1, sd), decreasing = TRUE))[1:n_top_regulons], ]
  
  col_fun <- colorRamp2(
    breaks = c(-2, 0, 2),  # Z-score 范围
    colors = c("white", "gray90", "black") 
  )
  # 绘制热图
  acell_scaled_data <- t(scale(t(regulonAUC)))
  pdf('G_SCENIC_Allcell_Heatmap.pdf',width = 5,height = 5)
  pheatmap(acell_scaled_data, 
           scale = "row", 
           clustering_method = "ward.D2",
           color = col_fun,
           show_colnames = FALSE,
           fontsize_row = 6)
  dev.off()
  
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_types <- G_sce$Celltype # 假设seurat对象中有celltype信息
  
  regulon_activity_by_celltype <- sapply(split(colnames(regulonAUC), cell_types), 
                                         function(cells) rowMeans(regulonAUC[,cells]))
  
  # 选择最具细胞类型特异性的regulon
  top_regulon_per_celltype <- apply(regulon_activity_by_celltype, 2, 
                                    function(x) names(x)[which.max(x)])
  cluster_to_celltype <- c(
    "Basophil" = "Basophil",
    "Eosinophil" = "Eosinophil",
    "Neutrophil" = "Neutrophil"
  )
  ###"#dc8e97","#e3d1db","#74a893","#e98741","#5ac6e9","#ebce8e","#e5c06e","#7587b1","#c7deef","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac","#64ae79"
  celltype_colors<- c(
    "Basophil" = "#dc8e97",
    "Eosinophil" = "#e3d1db",
    "Neutrophil" = "#5ac6e9"
  )

  heatmap_columns <- colnames(regulon_activity_by_celltype)
  
  column_ha <- HeatmapAnnotation(
    CellType = cluster_to_celltype[heatmap_columns],
    col = list(CellType = celltype_colors),  # 颜色映射
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  scaled_data <- t(scale(t(regulon_activity_by_celltype)))
  
  pdf('G_SCENIC_Celltype_Heatmap.pdf',width = 6,height = 7)
  Heatmap(scaled_data,
          name = "Regulon activity",
          top_annotation = column_ha,
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  # 计算每个regulon在不同细胞类型中的平均活性
  cell_Group <- G_sce$Group # 假设seurat对象中有celltype信息
  regulon_activity_by_Group <- sapply(split(colnames(regulonAUC), cell_Group), 
                                      function(cells) rowMeans(regulonAUC[,cells]))
  # 选择最具细胞类型特异性的regulon
  top_regulon_Group  <- apply(regulon_activity_by_Group, 2, 
                              function(x) names(x)[which.max(x)])
  #scaled_data_Group <- t(scale(t(regulon_activity_by_Group)))
  write.table(regulon_activity_by_Group,'G_SCENIC_Group.xls',sep = '\t') 
  pdf('G_SCENIC_Group_Heatmap.pdf',width = 6,height = 6)
  Heatmap(regulon_activity_by_Group,
          name = "Regulon activity",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 10))
  dev.off()
  write.table(L_sce@meta.data,'G_SCENIC_sce.xls',sep = '\t') 
  
  grn <- read.table('T_grn_output.tsv',sep = '\t',check.names = F,header = T) 
  pparg <- grn[which(grn$TF == 'PPARG'),]
  elf1 <- grn[which(grn$TF == 'ELF1'),]
  write.table(pparg,'TH1_PPARG_TF_Target.xls',sep = '\t')
  write.table(elf1,'TH1_ELF1_TF_Target.xls',sep = '\t')
  
  library(ggplot2)
  library(ggridges) # 用于绘制山峦图
  library(dplyr)    # 用于数据处理
  library(tidyr)    # 用于数据整理
  data <- read.table('G_SCENIC_sce.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  data <- data[,c(-1:-3,-5:-12)]
  data_long <- data %>%
    pivot_longer(cols = -Group, names_to = "TF", values_to = "AUC")
  
  ###NFE2L2(+) CEBPB(+)
  selected_TF <- "NFE2L2(+)"
  pdf('G_Ridges_Activity of NFE2L2(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
  
  selected_TF <- "CEBPB(+)"
  pdf('G_Ridges_Activity of CEBPB(+) Across Group.pdf',width = 5,height = 4)
  data_long %>% 
    filter(TF == selected_TF) %>%
    ggplot(aes(x = AUC, y = Group, fill = Group)) +
    geom_density_ridges(scale = 1.5, alpha = 0.7) +
    theme_ridges() +
    labs(title = paste("Activity of", selected_TF, "Across Group"),
         x = "AUCell score",
         y = "Group") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  dev.off()
}

###singlecell celltype heatmap
if(T){
  ###Ecell
  library(tidyverse)
  library(Seurat)
  library(zellkonverter)
  library(ClusterGVis)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(Seurat)
  scedata <- readRDS('Enormalize.sce_v4.rds')
  table(scedata@meta.data$Celltype2)
  scedata@active.ident <- as.factor(scedata@meta.data$Celltype2)
  Idents(scedata) <- scedata@meta.data[["Celltype2"]]
  pbmc.markers.all <- Seurat::FindAllMarkers(scedata, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
  write.table(pbmc.markers.all,'Ecell.markers.all.xls',sep = '\t')
  pbmc.markers.all <- read.table('Ecell.markers.all.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  pbmc.markers <- pbmc.markers.all %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
  st.data <- prepareDataFromscRNA(object = scedata,diffData = pbmc.markers,showAverage = TRUE)
  
  enrich <- enrichCluster(object = st.data,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.5,topn = 5, seed = 5201314)
  markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = F)]
  #visCluster(object = st.data, plot.type = "line")
  
  pdf('Ecell_scHeatmap.pdf',height = 12,width = 16,onefile = F)
  visCluster(object = st.data,
             plot.type = "both",
             column_names_rot = 70,
             show_row_dend = F,
             markGenes = markGenes,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             cluster.order = c(1:6),
             go.col = rep(jjAnno::useMyCol("stallion",n = 6),each = 5),
             add.bar = T)
  dev.off()
  
  ###Lcell
  library(tidyverse)
  library(Seurat)
  library(zellkonverter)
  library(ClusterGVis)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(Seurat)
  scedata <- readRDS('Lnormalize.sce_v4.rds')
  table(scedata@meta.data$Celltype)
  scedata@active.ident <- as.factor(scedata@meta.data$Celltype)
  Idents(scedata) <- scedata@meta.data[["Celltype"]]
  pbmc.markers.all <- Seurat::FindAllMarkers(scedata, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
  write.table(pbmc.markers.all,'Lcell.markers.all.xls',sep = '\t')
  pbmc.markers.all <- read.table('Lcell.markers.all.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  pbmc.markers <- pbmc.markers.all %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
  st.data <- prepareDataFromscRNA(object = scedata,diffData = pbmc.markers,showAverage = TRUE)
  
  enrich <- enrichCluster(object = st.data,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.5,topn = 5, seed = 5201314)
  markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = F)]
  #visCluster(object = st.data, plot.type = "line")
  
  pdf('Lcell_scHeatmap.pdf',height = 12,width = 16,onefile = F)
  visCluster(object = st.data,
             plot.type = "both",
             column_names_rot = 70,
             show_row_dend = F,
             markGenes = markGenes,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             cluster.order = c(1:10),
             go.col = rep(jjAnno::useMyCol("stallion",n = 10),each = 5),
             add.bar = T)
  dev.off()
  
  ###Mcell
  library(tidyverse)
  library(Seurat)
  library(zellkonverter)
  library(ClusterGVis)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(Seurat)
  scedata <- readRDS('Mnormalize.sce_v4.rds')
  table(scedata@meta.data$Celltype)
  scedata@active.ident <- as.factor(scedata@meta.data$Celltype)
  Idents(scedata) <- scedata@meta.data[["Celltype"]]
  pbmc.markers.all <- Seurat::FindAllMarkers(scedata, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
  write.table(pbmc.markers.all,'Mcell.markers.all.xls',sep = '\t')
  pbmc.markers.all <- read.table('Mcell.markers.all.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  pbmc.markers <- pbmc.markers.all %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
  st.data <- prepareDataFromscRNA(object = scedata,diffData = pbmc.markers,showAverage = TRUE)
  
  enrich <- enrichCluster(object = st.data,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.5,topn = 5, seed = 5201314)
  markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = F)]
  #visCluster(object = st.data, plot.type = "line")
  
  pdf('Mcell_scHeatmap.pdf',height = 12,width = 16,onefile = F)
  visCluster(object = st.data,
             plot.type = "both",
             column_names_rot = 70,
             show_row_dend = F,
             markGenes = markGenes,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             cluster.order = c(1:14),
             go.col = rep(jjAnno::useMyCol("stallion",n = 14),each = 5),
             add.bar = T)
  dev.off()

  ###Gcell
  library(tidyverse)
  library(Seurat)
  library(zellkonverter)
  library(ClusterGVis)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(Seurat)
  scedata <- readRDS('Nnormalize.sce_v4.rds')
  table(scedata@meta.data$Celltype)
  scedata@active.ident <- as.factor(scedata@meta.data$Celltype)
  Idents(scedata) <- scedata@meta.data[["Celltype"]]
  pbmc.markers.all <- Seurat::FindAllMarkers(scedata, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
  write.table(pbmc.markers.all,'Gcell.markers.all.xls',sep = '\t')
  pbmc.markers.all <- read.table('Gcell.markers.all.xls',sep = '\t',header = T,row.names = 1,check.names = F)
  pbmc.markers <- pbmc.markers.all %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
  st.data <- prepareDataFromscRNA(object = scedata,diffData = pbmc.markers,showAverage = TRUE)
  
  enrich <- enrichCluster(object = st.data,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.5,topn = 5, seed = 5201314)
  markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = F)]
  #visCluster(object = st.data, plot.type = "line")
  
  pdf('Gcell_scHeatmap.pdf',height = 12,width = 16,onefile = F)
  visCluster(object = st.data,
             plot.type = "both",
             column_names_rot = 70,
             show_row_dend = F,
             markGenes = markGenes,
             markGenes.side = "left",
             annoTerm.data = enrich,
             line.side = "left",
             cluster.order = c(1:3),
             go.col = rep(jjAnno::useMyCol("stallion",n = 3),each = 5),
             add.bar = T)
  dev.off()
}

###singleCell OR
if(T){
  do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                            meta.cluster = cellInfo.tb$meta.cluster,
                            colname.patient = "patient",
                            loc = cellInfo.tb$loc,
                            out.prefix,
                            pdf.width=3,
                            pdf.height=5,
                            verbose=0){
    ##input data 
    library(data.table)
    dir.create(dirname(out.prefix),F,T)
    
    cellInfo.tb = data.table(cellInfo.tb)
    cellInfo.tb$meta.cluster = as.character(meta.cluster)
    
    if(is.factor(loc)){
      cellInfo.tb$loc = loc
    }else{cellInfo.tb$loc = as.factor(loc)}
    
    loc.avai.vec <- levels(cellInfo.tb[["loc"]])
    count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
    freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)
    
    {
      count.dist.melt.ext.tb <- test.dist.table(count.dist)
      p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
      OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
      OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
      rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
    }
    
    sscVis::plotMatrix.simple(OR.dist.mtx,
                              out.prefix=sprintf("%s.OR.dist",out.prefix),
                              show.number=F,
                              waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                              exp.name=expression(italic(OR)),
                              z.hi=4,
                              palatte=viridis::viridis(7),
                              pdf.width = 4, pdf.height = pdf.height)
    if(verbose==1){
      return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                  "p.dist.tb"=p.dist.tb,
                  "OR.dist.tb"=OR.dist.tb,
                  "OR.dist.mtx"=OR.dist.mtx))
    }else{
      return(OR.dist.mtx)
    }
  }
  
  test.dist.table <- function(count.dist,min.rowSum=0)
  {
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
      this.row <- count.dist.melt.tb$rid[i]
      this.col <- count.dist.melt.tb$cid[i]
      this.c <- count.dist.melt.tb$count[i]
      other.col.c <- sum.col[this.col]-this.c
      this.m <- matrix(c(this.c,
                         sum.row[this.row]-this.c,
                         other.col.c,
                         sum(sum.col)-sum.row[this.row]-other.col.c),
                       ncol=2)
      res.test <- fisher.test(this.m)
      data.frame(rid=this.row,
                 cid=this.col,
                 p.value=res.test$p.value,
                 OR=res.test$estimate)
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    return(count.dist.melt.ext.tb)
  }

library(sscVis)
library(data.table)
library(grid)
library(cowplot)
library(ggrepel)
library(readr)
library(plyr)
library(ggpubr)
library(tidyverse)
library(viridis)
library(Seurat)
library(pheatmap)

  meta <- read.table('adata_E.obs.csv',sep = ',',header = T,check.names = F)
  meta$loc <- meta$Group
  meta$meta.cluster <- meta$Celltype2
  
  out.prefix <- "./Fig_OR"
  OR.immune.list <- do.tissueDist(cellInfo.tb=meta,
                                  out.prefix=sprintf("%s.Immune_cell",out.prefix),
                                  pdf.width=4,pdf.height=8,verbose=1
  )
  OR.immune.list
  
  a=OR.immune.list[["OR.dist.tb"]]
  a <- as.data.frame(a)
  rownames(a) <- a$rid
  a <- a[,-1]
  a <- na.omit(a)
  a
  
  b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
  b <- spread(b,key = "cid", value = "adj.p.value")
  b <- data.frame(b[,-1],row.names = b$rid)
  b <- b[rownames(a),]
  b
  
  col <- viridis(11,option = "D")
  b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
             ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                    ifelse(b<0.001&(a>1.5|a<0.5),"***",
                           ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                  ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))
  
  bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))
  
  pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
           display_numbers = b,number_color="black",fontsize_number=10,
           cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
           color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                     colorRampPalette(colors = col[6:11])(length(bk)/2)))
  
  meta <- read.table('adata_M.obs.csv',sep = ',',header = T,check.names = F)
  meta$loc <- meta$Group
  meta$meta.cluster <- meta$Celltype
  
  out.prefix <- "./Fig_OR"
  OR.immune.list <- do.tissueDist(cellInfo.tb=meta,
                                  out.prefix=sprintf("%s.Immune_cell",out.prefix),
                                  pdf.width=4,pdf.height=8,verbose=1
  )
  OR.immune.list
  
  a=OR.immune.list[["OR.dist.tb"]]
  a <- as.data.frame(a)
  rownames(a) <- a$rid
  a <- a[,-1]
  a <- na.omit(a)
  a
  
  b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
  b <- spread(b,key = "cid", value = "adj.p.value")
  b <- data.frame(b[,-1],row.names = b$rid)
  b <- b[rownames(a),]
  b
  
  col <- viridis(11,option = "D")
  b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
             ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                    ifelse(b<0.001&(a>1.5|a<0.5),"***",
                           ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                  ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))
  
  bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))
  
  pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
           display_numbers = b,number_color="black",fontsize_number=10,
           cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
           color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                     colorRampPalette(colors = col[6:11])(length(bk)/2)))

  meta <- read.table('adata_N.obs.csv',sep = ',',header = T,check.names = F)
  meta$loc <- meta$Group
  meta$meta.cluster <- meta$Celltype
  
  out.prefix <- "./Fig_OR"
  OR.immune.list <- do.tissueDist(cellInfo.tb=meta,
                                  out.prefix=sprintf("%s.Immune_cell",out.prefix),
                                  pdf.width=4,pdf.height=8,verbose=1
  )
  OR.immune.list
  
  a=OR.immune.list[["OR.dist.tb"]]
  a <- as.data.frame(a)
  rownames(a) <- a$rid
  a <- a[,-1]
  a <- na.omit(a)
  a
  
  b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
  b <- spread(b,key = "cid", value = "adj.p.value")
  b <- data.frame(b[,-1],row.names = b$rid)
  b <- b[rownames(a),]
  b
  
  col <- viridis(11,option = "D")
  b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
             ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                    ifelse(b<0.001&(a>1.5|a<0.5),"***",
                           ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                  ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))
  
  bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))
  
  pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
           display_numbers = b,number_color="black",fontsize_number=10,
           cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
           color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                     colorRampPalette(colors = col[6:11])(length(bk)/2)))
  
  
  }

###Cellchat
if(T){
  library(tidyverse)
  library(Seurat)
  library(CellChat)
  library(ComplexHeatmap)
  sce <- read.table('raw_matrix.csv.gz',sep = '\t',header = T,row.names = 1,check.names = F)
  metadata <- read.table('adata.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
  scedata <- CreateSeuratObject(counts = sce,meta.data = metadata,project = 'COPD')
  saveRDS(scedata,'scedata.rds')
  table(scedata$Group)
  sceCtrl <- subset(x = scedata,subset = Group == "Ctrl")
  data.input <- sceCtrl@assays$RNA$counts
  meta <- sceCtrl@meta.data
  cellchatCtrl <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
  cellchat <- cellchatCtrl
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "Cellchat_Ctrl_cell-cell_communications.all.csv")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  cellchatCtrl <- cellchat
  saveRDS(cellchatCtrl,'cellchat_Ctrl.rds')
  
  sceCOPD <- subset(x = scedata,subset = Group == "COPD")
  data.input <- sceCOPD@assays$RNA$counts
  meta <- sceCOPD@meta.data
  cellchatCOPD <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
  cellchat <- cellchatCOPD
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "Cellchat_COPD_cell-cell_communications.all.csv")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  cellchatCOPD <- cellchat
  saveRDS(cellchatCOPD,'cellchat_COPD.rds')
  
  cellchat.list <- list(Ctrl = cellchatCtrl,COPD = cellchatCOPD)
  cellchat <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = TRUE)
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  p <- gg1 + gg2
  p
  ggsave("Cellchat_Overview number strength.pdf",p,width = 6, height = 4)
  
  
  # 获取最大权重（用于统一比例）
  weight.max <- getMaxWeight(cellchat.list, attribute = c("idents", "count"))
  # 设置输出目录（可选）
  output_dir <- "CellChat_plots"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # 循环导出每个对象的 PDF
  for (i in 1:length(cellchat.list)) {
    # 设置 PDF 文件名（按组名命名）
    pdf_file <- file.path(output_dir, paste0("CellChat_Network_", names(cellchat.list)[i], ".pdf"))
    
    # 开启 PDF 设备
    pdf(pdf_file, width = 8, height = 6)  # 可调整宽高
    
    # 绘制网络图
    netVisual_circle(
      cellchat.list[[i]]@net$count,
      weight.scale = TRUE,
      label.edge = FALSE,
      edge.weight.max = weight.max[2],
      edge.width.max = 6,
      title.name = paste0("Number of interactions ~ ", names(cellchat.list)[i])
    )
    
    # 关闭 PDF 设备
    dev.off()
    
    # 提示保存成功
    message("Saved: ", pdf_file)
  }
  
  pathway.union <- union(cellchat.list[[1]]@netP$pathways,cellchat.list[[2]]@netP$pathways)
  
  pathway_scores_list <- list()
  
  for(i in 1:length(cellchat.list)) {
    # 提取通路活性得分矩阵
    score_mat <- cellchat.list[[i]]@netP$prob  # 或者使用其他适当的矩阵
    # 计算每个通路的平均得分
    pathway_scores_list[[i]] <- rowMeans(score_mat)
  }
  
  # 创建数据框计算变化
  pathway_changes <- data.frame(
    pathway = pathway.union
  )
  
  # 填充每个样本的得分
  for(i in 1:length(cellchat.list)) {
    pathway_changes[[paste0("sample", i)]] <- 
      pathway_scores_list[[i]][match(pathway.union, names(pathway_scores_list[[i]]))]
  }
  
  # 计算方差（变化程度）
  pathway_changes$variance <- apply(pathway_changes[, -1], 1, var, na.rm = TRUE)
  
  # 按方差排序，选择前40个
  top40_pathways <- pathway_changes %>%
    arrange(desc(variance)) %>%
    head(40) %>%
    pull(pathway)
  
  # 使用筛选后的通路绘制热图
  ht1 = netAnalysis_signalingRole_heatmap(cellchat.list[[1]], pattern = "all", signaling = top40_pathways,
                                          title = names(cellchat.list)[1], width = 6, height = 16)
  ht2 = netAnalysis_signalingRole_heatmap(cellchat.list[[2]], pattern = "all", signaling = top40_pathways,
                                          title = names(cellchat.list)[2], width = 6, height = 16)
  
  ht1 + ht2
  
  pdf("CellChat_plots/Cellchat_Ctrl_Overall signaling patterns.pdf", width = 8, height = 8)
  draw(ht1)  # 绘制 Heatmap
  dev.off()  # 关闭图形设备
  
  pdf("CellChat_plots/Cellchat_COPD_Overall signaling patterns.pdf", width = 6, height = 8)
  draw(ht2)  # 绘制 Heatmap
  dev.off()  # 关闭图形设备
  
  
  if(T){
    # 2. 设置图形参数以获得更好的视觉效果
    par(mfrow = c(1, 2), mar = c(2, 2, 4, 2))  # 调整边距
    
    # 3. 计算所有组的最大交互数以统一比例
    max_count <- max(sapply(cellchat.list, function(x) max(x@net$count)))
    
    # 4. 绘制三个组的网络图
    for (i in 1:length(cellchat.list)) {
      # 获取当前组的细胞类型数量
      groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
      
      # 绘制网络图
      netVisual_circle(cellchat.list[[i]]@net$count, 
                       vertex.weight = groupSize,
                       weight.scale = TRUE, 
                       label.edge = FALSE,
                       edge.weight.max = max_count,  # 使用统一的最大值
                       title.name = paste0(names(cellchat.list)[i], 
                                           "\n(Total: ", sum(cellchat.list[[i]]@net$count), ")"))
    }
    # 5. 重置图形参数
    par(mfrow = c(1, 1))
  }
  
  if(T){
    # 2. 设置图形参数以获得更好的视觉效果
    par(mfrow = c(1, 2), mar = c(2, 2, 4, 2))  # 调整边距
    
    # 3. 计算所有组的最大交互数以统一比例
    max_count <- max(sapply(cellchat.list, function(x) max(x@net$count)))
    
    # 4. 绘制三个组的网络图
    for (i in 1:length(cellchat.list)) {
      # 获取当前组的细胞类型数量
      groupSize <- as.numeric(table(cellchat.list[[i]]@idents))
      
      # 绘制网络图
      netVisual_circle(cellchat.list[[i]]@net$weight, 
                       vertex.weight = groupSize,
                       weight.scale = TRUE, 
                       label.edge = FALSE,
                       edge.weight.max = max_count,  # 使用统一的最大值
                       title.name = paste0(names(cellchat.list)[i], 
                                           "\n(Total: ", sum(cellchat.list[[i]]@net$weight), ")"))
    }
    # 5. 重置图形参数
    par(mfrow = c(1, 1))
  }
  
  ###Ctrl.vs.trial
  if(T){
    cellchat.list <- list(Ctrl = cellchatCtrl,COPD = cellchatCOPD)
    cellchat <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = TRUE)   
    par(mfrow = c(1,2))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    p1 <- h1 + h2
    pdf("CellChat_plots/Cellchat_Ctrl.vs.COPD_netVisual_heatmap.pdf", width = 10, height = 8)
    draw(p1)  # 绘制 Heatmap
    dev.off()  # 关闭图形设备
    
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T,do.stat = TRUE)
    ggsave("CellChat_plots/Cellchat_Ctrl.vs.trial.stacked.Compare_pathway_strengh.pdf",gg1,width = 6, height = 8)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F,do.stat = TRUE)
    ggsave("CellChat_plots/Cellchat_Ctrl.vs.trial.Barplot.Compare_pathway_strengh.pdf",gg2,width = 6, height = 8)
  }
}

if(T){
  metadata <- read.table('adata.obs.csv',sep = ',',header = T,check.names = F)
  head(metadata)
  library(dplyr)
  # 计算每个样本中各种细胞类型的百分比
  celltype_percent <- metadata %>%
    group_by(sample, Celltype) %>%
    summarise(count = n()) %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # 查看结果
  print(celltype_percent, n = Inf)
  
  library(tidyr)
  df <- celltype_percent[,-3]
  wide_df <- df %>%
    pivot_wider(
      names_from = Celltype,   # 将 Celltype 列中的值作为新列名
      values_from = percentage # 将 percentage 列的值填入对应单元格
    )
  
  # 查看转换后的宽表格
  head(wide_df)
  
  write.table(wide_df,'Allcelltype_samle_percent.xls',sep = '\t')
}

if(T){
  metadata <- read.table('adata_E.obs.csv',sep = ',',header = T,check.names = F)
  head(metadata)
  library(dplyr)
  # 计算每个样本中各种细胞类型的百分比
  celltype_percent <- metadata %>%
    group_by(sample, Celltype) %>%
    summarise(count = n()) %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # 查看结果
  print(celltype_percent, n = Inf)
  
  library(tidyr)
  df <- celltype_percent[,-3]
  wide_df <- df %>%
    pivot_wider(
      names_from = Celltype,   # 将 Celltype 列中的值作为新列名
      values_from = percentage # 将 percentage 列的值填入对应单元格
    )
  
  # 查看转换后的宽表格
  head(wide_df)
  
  write.table(wide_df,'Ecelltype_samle_percent.xls',sep = '\t')
}

if(T){
  metadata <- read.table('adata_L.obs.csv',sep = ',',header = T,check.names = F)
  head(metadata)
  library(dplyr)
  # 计算每个样本中各种细胞类型的百分比
  celltype_percent <- metadata %>%
    group_by(sample, Celltype) %>%
    summarise(count = n()) %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # 查看结果
  print(celltype_percent, n = Inf)
  
  library(tidyr)
  df <- celltype_percent[,-3]
  wide_df <- df %>%
    pivot_wider(
      names_from = Celltype,   # 将 Celltype 列中的值作为新列名
      values_from = percentage # 将 percentage 列的值填入对应单元格
    )
  
  # 查看转换后的宽表格
  head(wide_df)
  
  write.table(wide_df,'Lcelltype_samle_percent.xls',sep = '\t')
}

if(T){
  metadata <- read.table('adata_M.obs.csv',sep = ',',header = T,check.names = F)
  head(metadata)
  library(dplyr)
  # 计算每个样本中各种细胞类型的百分比
  celltype_percent <- metadata %>%
    group_by(sample, Celltype) %>%
    summarise(count = n()) %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # 查看结果
  print(celltype_percent, n = Inf)
  
  library(tidyr)
  df <- celltype_percent[,-3]
  wide_df <- df %>%
    pivot_wider(
      names_from = Celltype,   # 将 Celltype 列中的值作为新列名
      values_from = percentage # 将 percentage 列的值填入对应单元格
    )
  
  # 查看转换后的宽表格
  head(wide_df)
  
  write.table(wide_df,'Mcelltype_samle_percent.xls',sep = '\t')
}

if(T){
  metadata <- read.table('adata_N.obs.csv',sep = ',',header = T,check.names = F)
  head(metadata)
  library(dplyr)
  # 计算每个样本中各种细胞类型的百分比
  celltype_percent <- metadata %>%
    group_by(sample, Celltype) %>%
    summarise(count = n()) %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # 查看结果
  print(celltype_percent, n = Inf)
  
  library(tidyr)
  df <- celltype_percent[,-3]
  wide_df <- df %>%
    pivot_wider(
      names_from = Celltype,   # 将 Celltype 列中的值作为新列名
      values_from = percentage # 将 percentage 列的值填入对应单元格
    )
  
  # 查看转换后的宽表格
  head(wide_df)
  
  write.table(wide_df,'Gcelltype_samle_percent.xls',sep = '\t')
}

if(T){
  library(Seurat)
  df <- read.table('新建文件夹/GSM5575235_NGS19_K360_Dim_LBA_Hum12.counts_tbl.csv.gz',sep = ',',header = T,row.names = 1,check.names = F)
  scedata <- CreateSeuratObject(counts = df)
  # 加载必要包
  library(Seurat)
  library(Matrix)
  
  counts_matrix <- GetAssayData(scedata, assay = "RNA", layer = "counts")
  
  sample_ids <- scedata$orig.ident
  gene_names <- rownames(counts_matrix)
  gene_ids <- gene_names  # 若无可用的外部ID，则与基因名相同
  feature_type <- rep("Gene Expression", length(gene_names))
  features_df <- data.frame(GeneID = gene_ids,
                            GeneName = gene_names,
                            FeatureType = feature_type)
  
  # 为每个样本导出三文件到单独文件夹（建议每个样本一个子目录）
  unique_samples <- unique(sample_ids)
  
  for (samp in unique_samples) {
    # 创建样本专属文件夹
    dir_name <- samp
    if (!dir.exists(dir_name)) dir.create(dir_name)
    
    # 选择当前样本的细胞
    cells_use <- which(sample_ids == samp)
    sub_counts <- counts_matrix[, cells_use, drop = FALSE]  # 稀疏矩阵子集
    
    # 准备 barcodes 向量（细胞条形码）
    barcodes <- colnames(sub_counts)  # 通常是原始细胞名，如 "AAACCCATCTACAGGT-1_1"
    
    # 导出 matrix.mtx（需注意 10x 格式是 行列转置的？不，10x 的标准是基因×细胞，列索引为细胞）
    # Matrix 包的 writeMM 写入的是 Matrix Market 格式，列序与矩阵一致，直接使用即可。
    # 注意：10x 的 matrix.mtx 文件头三行注释，且行列数量需正确。
    mtx_file <- file.path(dir_name, "matrix.mtx")
    writeMM(sub_counts, file = mtx_file)
    # 修正文件头：Matrix::writeMM 输出的头为标准形式，但10x 要求第一行为 %%MatrixMarket matrix coordinate real general
    # writeMM 已生成正确头，无需额外修改。
    
    # 导出 features.tsv
    features_file <- file.path(dir_name, "features.tsv")
    write.table(features_df, file = features_file, sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # 导出 barcodes.tsv
    barcodes_file <- file.path(dir_name, "barcodes.tsv")
    writeLines(barcodes, con = barcodes_file)
    
    cat("已完成样本", samp, " : ", ncol(sub_counts), "细胞,", nrow(sub_counts), "基因\n")
  }
}