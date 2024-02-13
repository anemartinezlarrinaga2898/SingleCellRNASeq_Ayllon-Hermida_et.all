# SCRIPT: Seurat Basic Pipeline 

# FECHA: 03/03/2022
# AUTOR: ANE MARTINEZ LARRINAGA 

# PORTILLO SAMPLE A0059

################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/2_CELLULAR_SYSTEMS_GENOMICS/5_PORTILLO")

library(Seurat)
library(tidyverse)
library(DoubletFinder)

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

path.guardar <- "0.1_SingleCellAnalysis"
path.files <- "0.0_CellRanger"
file.name <- "A0059"

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

ReadFiles <- function(path,file.name){
  path <- paste(path,file.name,"filtered_feature_bc_matrix",sep="/")
  sample <- Read10X(paste(path.files,file.name,"filtered_feature_bc_matrix",sep="/"))
  return(sample) 
}

CreateSeuratObejct <- function(matrix,nombre.projecto,condicion.estudio){
  data <-CreateSeuratObject(counts = matrix, project = nombre.projecto, min.cells = 3)
  data@meta.data[["ID"]] <- nombre.projecto
  data@meta.data[["Sample.Type"]] <- condicion.estudio
  return(data)
}

QualityControlCheck <- function(data,pattern.mitocondrial,pattern.ribosomal){
  data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
  data$percent.mt <- PercentageFeatureSet(data, pattern = pattern.mitocondrial)
  data$percent.rb <- PercentageFeatureSet(data, pattern = pattern.ribosomal)
  
  metadata <- data@meta.data  %>%
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  
  data@meta.data <- metadata
  return(data)
}

GraficosQualityControl <- function(data){
  metadata <- data@meta.data
  # Visualize the number of cell counts per sample
  
  NumberCellsPerSample <-metadata %>% 
    ggplot(aes(x=Sample.Type, fill=ID)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  
  # Visualize the number UMIs/transcripts per cell
  
  NumberOfUMIPerCell <- metadata %>% 
    ggplot(aes(color=Sample.Type, x=nUMI, fill= Sample.Type)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("log10 Cell density")+
    ggtitle("NumberOfUMIPerCell")
  
  # Visualize the distribution of genes detected per cell via histogram
  
  GenesDetectedPerCell <- metadata %>% 
    ggplot(aes(color=Sample.Type, x=nGene, fill= Sample.Type)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10()+
    ggtitle("GenesDetectedPerCell")
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  
  Complexity <- metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = Sample.Type, fill=Sample.Type)) +
    geom_density(alpha = 0.2) +
    theme_classic()+
    ggtitle("Complexity")
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  
  MitoPerCell <- metadata %>% 
    ggplot(aes(color=Sample.Type, x=percent.mt, fill=Sample.Type)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic()+
    ggtitle("MitoPerCell")
  
  RiboPerCell <- metadata %>% 
    ggplot(aes(color=Sample.Type, x=percent.rb, fill=Sample.Type)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic()+
    ggtitle("MitoPerCell")
  
  # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  
  Correlation <- metadata %>% 
    ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    ggtitle("Correlation")+
    theme_classic()+
    facet_wrap(~Sample.Type) 
  
  ViolinPlots.logscale <- VlnPlot(data, features = c("nGene", "nUMI", "percent.mt","percent.rb"), ncol = 4, log = T)
  ViolinPlots.No.logscale <- VlnPlot(data, features = c("nGene", "nUMI", "percent.mt","percent.rb"), ncol = 4, log = F)
  
  lista.graficos <- list(NumberCellsPerSample=NumberCellsPerSample,
                         NumberOfUMIPerCell=NumberOfUMIPerCell,
                         GenesDetectedPerCell=GenesDetectedPerCell,
                         Complexity=Complexity,
                         MitoPerCell=MitoPerCell,
                         RiboPerCell=RiboPerCell,
                         Correlation=Correlation,
                         ViolinPlots.logscale=ViolinPlots.logscale,
                         ViolinPlots.No.logscale=ViolinPlots.No.logscale)
  return(lista.graficos)
}

EliminacionDoublets <- function(data){
  require(DoubletFinder)
  # pre-process standard workflow.................................................
  data <- NormalizeData(object = data)
  data <- FindVariableFeatures(object = data)
  data <- ScaleData(object = data)
  data <- RunPCA(object = data)
  
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  dim.final <- min(co1, co2)
  
  ElbowPlot <- ElbowPlot(data)
  
  data <- FindNeighbors(object = data, dims = 1:dim.final)
  data <- FindClusters(object = data)
  data <- RunUMAP(object = data, dims = 1:dim.final)
  
  ## pK Identification (no ground-truth)..........................................
  # No ground truth means that we dont know if a cell is a single o a doublet
  # There are experimental methods to know this, but normally people dont know this
  
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.final, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  
  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  ## Homotypic Doublet Proportion Estimate........................................
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.076*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # run doubletFinder ............................................................
  data <- doubletFinder_v3(data,
                           PCs = 1:dim.final,
                           pN = 0.25,
                           pK = pK,
                           nExp = nExp_poi.adj,
                           reuse.pANN = FALSE,
                           sct = FALSE)
  
  idx.colnames.doubletValue <- stringr::str_detect(string = colnames(data@meta.data),pattern = "^pANN_")
  idx.colnames.doubletClasification <- stringr::str_detect(string = colnames(data@meta.data),pattern = "^DF.classifications")
  
  colnames(data@meta.data)[idx.colnames.doubletValue] <- "DoubletValue"
  colnames(data@meta.data)[idx.colnames.doubletClasification] <- "DoubletClassification"
  
  # visualize doublets............................................................
  Lista.Vuelta <- list(data=data,
                       ElbowPlot=ElbowPlot,
                       dim.final=dim.final)
  return(Lista.Vuelta)
}


FuncionColores <- function(n){
  require(RColorBrewer)
  n <- n
  qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
}
col <- FuncionColores(15)

Calculo.Dimensiones.PCA <- function(data){
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  dim.final <- min(co1, co2)
  return(dim.final)
}

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
################################################################################

# 1) Create Seurat Object ......................................................

sample <- ReadFiles(path=path.files,file.name = file.name)

seurat.object <- CreateSeuratObejct(matrix = sample,
                                    nombre.projecto=file.name,
                                    condicion.estudio="A0059")

# 2) Quality Control Check: ....................................................
# Calculate quality metrix in order to filter afterwards
seurat.object <- QualityControlCheck(seurat.object,pattern.mitocondrial="^mt-",pattern.ribosomal="^rb[sl]")
# Now you estimate the graphs. 
graphs.quality.control <- GraficosQualityControl(data=seurat.object)

pdf(paste(path.guardar,"A0059_Graficos.QualityControl.pdf",sep="/"),width = 15,height = 20)
graphs.quality.control
dev.off()

# 3) Filter the samples ........................................................

data.filtered <- subset(seurat.object, subset = nGene > 300 & log10GenesPerUMI > 0.8)
data.filtered <- subset(data.filtered, subset = percent.mt < 10)

ViolinPlots.logscale <- VlnPlot(data.filtered, features = c("nGene", "nUMI", "percent.mt","percent.rb"), ncol = 4, log = T)
ViolinPlots.No.logscale <- VlnPlot(data.filtered, features = c("nGene", "nUMI", "percent.mt","percent.rb"), ncol = 4, log = F)

graphs.quality.control.after.filtering <- list(ViolinPlots.logscale,ViolinPlots.No.logscale)

pdf(paste(path.guardar,"A0059_Graficos.After.QualityControl.pdf",sep="/"),width = 15,height = 20)
graphs.quality.control.after.filtering
dev.off()

# 4) Identify Doublets .........................................................

data.filtered <- EliminacionDoublets(data = data.filtered)
Umap.doublet <- DimPlot(data.filtered$data, 
                        reduction = "umap", 
                        group.by = "DoubletClassification", 
                        label = TRUE, 
                        label.size = 4)+ggtitle("DoubletIdentificacion")
ggsave(Umap.doublet,filename = paste(path.guardar,"A0059_DoubletIdentification.pdf",sep = "/"),width = 10,height = 10)

data.filtered <- data.filtered$data

proportions.doublets <- as.data.frame(table(data.filtered$DoubletClassification))
write_csv(proportions.doublets,paste(path.guardar,"A0059_Proportions.Doublets.csv",sep="/"))

# 5) Seurat Pipeline: ..........................................................

# ..............................................................................
data <- data.filtered
# ..............................................................................

data <- NormalizeData(object = data)

gene.info.distribution <- summary(Matrix::colSums(data@assays$RNA@counts[,]>0))
hvg.number <- round(gene.info.distribution[4]+100)
data <- FindVariableFeatures(object = data,selection.method = "vst", nfeatures = hvg.number)
data<- ScaleData(object = data)
data<- RunPCA(object = data)

ElbowPlot <- ElbowPlot(data,ndims = 30)
ggsave(ElbowPlot,filename = paste(path.guardar,"A0059_ElbowPlot.png",sep = "/"),width = 10,height = 10)

dimensiones.PCA <- Calculo.Dimensiones.PCA(data = data)
print(dimensiones.PCA)

data <- FindNeighbors(data, dims = 1:dimensiones.PCA)
data <- FindClusters(data, resolution = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5))

data <- RunUMAP(data, reduction.use = "pca", dims = 1:dimensiones.PCA)

# Draw UMAPs by different resolutions: -----------------------------------------

umap.res.0.1 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.0.1",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.0.3 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.0.3",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.0.5 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.0.5",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.0.7 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.0.7",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.0.9 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.0.9",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.1.1 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.1.1",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.1.3 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.1.3",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

umap.res.1.5 <- DimPlot(data,
                        reduction = "umap",
                        group.by = "RNA_snn_res.1.5",
                        label = TRUE, 
                        label.size = 8,
                        cols =col,pt.size = 1)+NoLegend()

Umaps.totales <- list(umap.res.0.1,umap.res.0.3,umap.res.0.5,umap.res.0.7,umap.res.0.9,umap.res.1.1,umap.res.1.3,umap.res.1.5)

pdf(paste(path.guardar,"A0059_Umaps.DifferentResolutions.pdf",sep="/"),width = 15,height = 20)
Umaps.totales
dev.off()

# Estimate the number of cluster in each resolution 

metadata <- data@meta.data
idx.colnames.resolution <- stringr::str_detect(string = colnames(data@meta.data),pattern = "^RNA_snn")

metadata.resolution <- metadata[,idx.colnames.resolution]

number.cluster.for.each.resolution <-vector(length = ncol(metadata.resolution))
names(number.cluster.for.each.resolution) <- c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5)

i <- 1
for(i in seq(1,ncol(metadata.resolution),1)){
  res <- metadata.resolution[,i]
  unique.clusters <- as.double(unique(res))
  number.cluster.for.each.resolution[i] <- max(unique.clusters)
}

number.cluster.for.each.resolution <- as.data.frame(number.cluster.for.each.resolution)
write_csv(number.cluster.for.each.resolution,paste(path.guardar,"A0059_NumberOfClustersInEachResolution.csv",sep="/"))

saveRDS(data,paste(path.guardar,"A0059_SeuratObjectProcessed.rds",sep="/"))

