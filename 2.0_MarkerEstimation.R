# SCRIPT; Con los request que te piden el portillo 

# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 17-03-2023

################################################################################

directory <- setwd("5_PORTILLO")

library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(foreach)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
col <-  getPalette(8)

################################################################################

data.fibros <- readRDS("0.1_SingleCellAnalysis/FibroblastData.rds")
data.bela <- readRDS("0.1_SingleCellAnalysis/BelaData.rds")

# Obtener los genes diferencialmente expresados entre clusters 

  # Fibroblast
DimPlot(data.fibros, group.by = "RNA_snn_res.0.3",cols = col,pt.size = 2)&NoAxes()&ggtitle("Fibroblast_Res_0.3")
data.fibros <- SetIdent(data.fibros,value = "RNA_snn_res.0.3")
Markers.Clusters.Fibros <- FindAllMarkers(data.fibros)
clusters <- unique(Markers.Clusters.Fibros$cluster)
clusters <- clusters[order(clusters,decreasing = F)]
clusters.names <- paste("Cluster",clusters,sep="_")

Ordenar.Genes <- foreach::foreach(i=1:length(clusters),.final =function(x)setNames(x,clusters.names))%do%{
  c <- clusters[i]
  idx.c <- which(Markers.Clusters.Fibros$cluster==c)
  markers.cluster <- Markers.Clusters.Fibros[idx.c,]
  markers.cluster <- markers.cluster[order(markers.cluster$avg_log2FC,decreasing = T),]
  markers.cluster <- markers.cluster[1:50,]
}

Markers.Clusters.Fibros <- do.call(rbind,Ordenar.Genes)
Markers.Clusters.Fibros <- na.omit(Markers.Clusters.Fibros)
rownames(Markers.Clusters.Fibros) <- NULL

  # Bela 
DimPlot(data.bela, group.by = "RNA_snn_res.0.3",cols = col,pt.size = 2)&NoAxes()&ggtitle("Fibroblast_Res_0.3")
data.bela <- SetIdent(data.bela,value = "RNA_snn_res.0.3")
Markers.Clusters.Bela <- FindAllMarkers(data.bela)
clusters <- unique(Markers.Clusters.Bela$cluster)
clusters <- clusters[order(clusters,decreasing = F)]
clusters.names <- paste("Cluster",clusters,sep="_")

Ordenar.Genes <- foreach::foreach(i=1:length(clusters),.final =function(x)setNames(x,clusters.names))%do%{
  c <- clusters[i]
  idx.c <- which(Markers.Clusters.Bela$cluster==c)
  markers.cluster <- Markers.Clusters.Bela[idx.c,]
  markers.cluster <- markers.cluster[order(markers.cluster$avg_log2FC,decreasing = T),]
  markers.cluster <- markers.cluster[1:50,]
}

Markers.Clusters.Bela <- do.call(rbind,Ordenar.Genes)
Markers.Clusters.Bela <- na.omit(Markers.Clusters.Bela)
rownames(Markers.Clusters.Bela) <- NULL

# Join the data of the cluster maresr

Lista.Markers.Cluster <- list(Markers.Clusters.Fibros,Markers.Clusters.Bela)
names(Lista.Markers.Cluster) <- c("Fibroblast","Bela")

# Save the data 
openxlsx::write.xlsx(Lista.Markers.Cluster, file = '0.0_Request/Markers_ByCluster.xlsx') 

# Obtener los markers entre tratamiento y no 

  # Fibros
data.fibros <- SetIdent(data.fibros,value = "Phenotype")
Markers.Fenotype.Fibros <- FindAllMarkers(data.fibros)

  # Bela
data.bela <- SetIdent(data.bela,value = "Phenotype")
Markers.Fenotype.Bela <- FindAllMarkers(data.bela)

Lista.Markers.Phenotype <- list(Markers.Fenotype.Fibros,Markers.Fenotype.Bela)
names(Lista.Markers.Phenotype) <- c("Fibroblast","Bela")

# Save the data 
openxlsx::write.xlsx(Lista.Markers.Phenotype, file = '0.0_Request/Markers_ByPhenotype.xlsx') 

# The plots with the genes of interest 

genes.fibros <- c("CYTH3","CDH4","CADM1","CD44","ITGA1","ICAM1","ITGA3","THBS1","VCAM1","FN1","AMIGO2","CD36")
genes.belas <- c("IFIT1","HBG1","HBD","HBA2","IFIT2","HBG2","IFIT3","NTS","GOS2","TIMP3","GATA2","GATA1","ARID3A","ALAS2","TAL1","ALAS1")

getPalette <-  colorRampPalette(brewer.pal(3, "RdYlBu"))
col <-  getPalette(3)
col <- viridis::viridis(2)

DotPlot.Fibros <- DotPlot(data.fibros,features =genes.fibros,cols = col)&
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20))&
  coord_flip()
ggsave(filename = "0.0_Request/DotPlotFibros.pdf",width = 10,height = 10,plot = DotPlot.Fibros)

DotPlot.Bela <- DotPlot(data.bela,features = genes.belas,cols = col)&
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20))&
  coord_flip()
ggsave(filename = "0.0_Request/DotPlot.Bela.pdf",width = 10,height = 10,plot = DotPlot.Bela)







