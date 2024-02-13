# SCRIPT: Generar los objetos correspondientes para las celulas tratadas y las no tratadas
# de ambos cultivos 

# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 07-03-2023

################################################################################
directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/2_CELLULAR_SYSTEMS_GENOMICS/5_PORTILLO")

library(Seurat)
library(tidyverse)
library(scGate)
library(RColorBrewer)
library(patchwork)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
cols <-  getPalette(7)

SeuratPipeline <- function(data){
  gene.info.distribution <- summary(Matrix::colSums(data@assays$RNA@counts[,]>0))
  hvg.number <- round(gene.info.distribution[4]+100)
  data <- FindVariableFeatures(object = data,selection.method = "vst", nfeatures = hvg.number)
  data<- ScaleData(object = data)
  data<- RunPCA(object = data)
  
  ElbowPlot <- ElbowPlot(data,ndims = 30)

  dimensiones.PCA <- Calculo.Dimensiones.PCA(data = data)
  print(dimensiones.PCA)
  
  data <- FindNeighbors(data, dims = 1:dimensiones.PCA)
  data <- FindClusters(data, resolution = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5))
  data <- RunUMAP(data, reduction.use = "pca", dims = 1:dimensiones.PCA)
}

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
################################################################################

data.59 <- readRDS("0.1_SingleCellAnalysis/A0059_SeuratObjectProcessed.rds")
data.60 <- readRDS("0.1_SingleCellAnalysis/A0060_SeuratObjectProcessed.rds")

# Add the metadata to know which one is treated and which one is not 
data.59$Phenotype <- "NotTreated"
data.60$Phenotype <- "Treated"

genes.totales.59 <- rownames(data.59@assays$RNA@counts)
genes.totales.60 <- rownames(data.60@assays$RNA@counts)

# Selection of the fibroblast markers 
genes.markers <- c("Pdgfra","Col3a1","Col1a2","Abca8a","Dpt","S100a10","Gsn","Nav1","Cd34","Entpd2","Lum","Dcn")
genes.markers <- toupper(genes.markers)

gene.intersection.59 <- intersect(genes.totales.59,genes.markers)
gene.intersection.60 <- intersect(genes.totales.60,genes.markers)

cells.id.59 <- rownames(data.59@meta.data)
cells.id.60 <- rownames(data.60@meta.data)

# Selection of the Sample 59 ...................................................
fibroblast_model <- gating_model(name = "Fibroblast", signature = gene.intersection.59)
selection.59 <- scGate(data = data.59, model = fibroblast_model)

selection.cells.59 <- DimPlot(selection.59, group.by = "is.pure",cols = c("#64AAD2","#D5D5D5"))&NoAxes()
ggsave(filename = "0.1_SingleCellAnalysis/DimPlot_CellSelection_59.pdf",width = 10,height = 10,plot = selection.cells.59)

cells.selected.59 <- which(selection.59$is.pure=="Pure")

ID.cells.selected.59 <- rownames(data.59@meta.data)[cells.selected.59]
not.selected.59 <- setdiff(cells.id.59,ID.cells.selected.59)

fibroblast.59 <- subset(data.59,cells=ID.cells.selected.59)
bela.59 <- subset(data.59,cells=not.selected.59)

# Selection of the Sample 60 ...................................................
fibroblast_model <- gating_model(name = "Fibroblast", signature = gene.intersection.60)
selection.60 <- scGate(data = data.60, model = fibroblast_model)

selection.cells.60 <- DimPlot(selection.60, group.by = "is.pure",cols = c("#64AAD2","#D5D5D5"))&NoAxes()
ggsave(filename = "0.1_SingleCellAnalysis/DimPlot_CellSelection_60.pdf",width = 10,height = 10,plot = selection.cells.60)

cells.selected.60 <- which(selection.60$is.pure=="Pure")

ID.cells.selected.60 <- rownames(data.60@meta.data)[cells.selected.60]
not.selected.60 <- setdiff(cells.id.60,ID.cells.selected.60)

fibroblast.60 <- subset(data.60,cells=ID.cells.selected.60)
bela.60 <- subset(data.60,cells=not.selected.60)

# Generar Objeto fibros

data.fibroblastos <- merge(x=fibroblast.59,y=fibroblast.60)
data.bela <- merge(x=bela.59,y=bela.60)

# Seurat Pipeline for the merge object 

data.fibroblastos <- SeuratPipeline(data.fibroblastos)
data.bela <- SeuratPipeline(data.bela)

data.fibroblastos <- SetIdent(data.fibroblastos,value="RNA_snn_res.0.3")
data.bela <- SetIdent(data.bela,value="RNA_snn_res.0.3")

PhenotypeFibrosPlot <- DimPlot(data.fibroblastos,group.by = "Phenotype",cols = alpha(c("#6495BF","#F4D166"),0.66),pt.size = 3)&NoAxes()&ggtitle("Phenotype")
PhenotypeBela <- DimPlot(data.bela,group.by = "Phenotype",cols = alpha(c("#6495BF","#F4D166"),0.66),pt.size = 3)&NoAxes()&ggtitle("Phenotype")

ClustersFibros <- DimPlot(data.fibroblastos,cols =cols,pt.size = 3)&NoAxes()&ggtitle("Clusters")
ClustersBela <- DimPlot(data.bela,cols =cols,pt.size = 3)&NoAxes()&ggtitle("Clusters")

Fibros.Plots <- wrap_plots(list(ClustersFibros,PhenotypeFibrosPlot), ncol=2, guides="collect")
Bela.Plots <- wrap_plots(list(ClustersBela,PhenotypeBela), ncol=2, guides="collect")

ggsave(filename = "0.1_SingleCellAnalysis/Fibros.Plots.pdf",width = 25,height = 10,plot = Fibros.Plots)
ggsave(filename = "0.1_SingleCellAnalysis/Bela.Plots.pdf",width = 25,height = 10,plot = Bela.Plots)

# Removing doublets 

DoubletsFibros <- DimPlot(data.fibroblastos,group.by = "DoubletClassification",cols = alpha(c("#515A66","#7BC16E"),0.66),pt.size = 3)&NoAxes()&ggtitle("Phenotype")
DoubletsBelas <- DimPlot(data.bela,group.by = "DoubletClassification",cols = alpha(c("#515A66","#7BC16E"),0.66),pt.size = 3)&NoAxes()&ggtitle("Phenotype")

DoubletsFibros <- wrap_plots(list(ClustersFibros,DoubletsFibros), ncol=2, guides="collect")
DoubletsBelas <- wrap_plots(list(ClustersBela,DoubletsBelas), ncol=2, guides="collect")

ggsave(filename = "0.1_SingleCellAnalysis/Doublets_Fibros.Plots.pdf",width = 25,height = 10,plot = DoubletsFibros)
ggsave(filename = "0.1_SingleCellAnalysis/Doublets_Bela.Plots.pdf",width = 25,height = 10,plot = DoubletsBelas)

data.fibroblastos <- subset(data.fibroblastos,subset=DoubletClassification=="Singlet")
data.bela <- subset(data.bela,subset=DoubletClassification=="Singlet")

data.fibroblastos <- SeuratPipeline(data.fibroblastos)
data.bela <- SeuratPipeline(data.bela)

data.fibroblastos <- SetIdent(data.fibroblastos,value="RNA_snn_res.0.3")
data.bela <- SetIdent(data.bela,value="RNA_snn_res.0.3")

PhenotypeFibrosPlot <- DimPlot(data.fibroblastos,group.by = "Phenotype",cols = alpha(c("#6495BF","#F4D166"),0.66),pt.size = 3)&NoAxes()&ggtitle("Phenotype")
PhenotypeBela <- DimPlot(data.bela,group.by = "Phenotype",cols = alpha(c("#6495BF","#F4D166"),0.66),pt.size = 3)&NoAxes()&ggtitle("Phenotype")

ClustersFibros <- DimPlot(data.fibroblastos,cols =cols,pt.size = 3)&NoAxes()&ggtitle("Clusters")
ClustersBela <- DimPlot(data.bela,cols =cols,pt.size = 3)&NoAxes()&ggtitle("Clusters")

Fibros.Plots <- wrap_plots(list(ClustersFibros,PhenotypeFibrosPlot), ncol=2, guides="collect")
Bela.Plots <- wrap_plots(list(ClustersBela,PhenotypeBela), ncol=2, guides="collect")

ggsave(filename = "0.1_SingleCellAnalysis/Fibros_NoDoublets.Plots.pdf",width = 25,height = 10,plot = Fibros.Plots)
ggsave(filename = "0.1_SingleCellAnalysis/Bela_NoDoublets.Plots.pdf",width = 25,height = 10,plot = Bela.Plots)

saveRDS(data.fibroblastos,"0.1_SingleCellAnalysis/FibroblastData.rds")
saveRDS(data.bela,"0.1_SingleCellAnalysis/BelaData.rds")
