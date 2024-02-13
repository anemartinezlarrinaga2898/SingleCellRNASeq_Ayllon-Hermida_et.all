# SingleCellRNASeq_Ayllon-Hermida_et.all

**Summary paper:** 

Recent studies indicate that the human spleen contains over 95% of the total parasite biomass during chronic asymptomatic infections caused by Plasmodium vivax. Previous studies have demonstrated that extracellular vesicles (EVs) secreted from infected reticulocytes facilitate binding to human spleen fibroblasts (hSFs) and identified parasite genes whose expression was dependent on an intact spleen. Here, we characterize the P. vivax spleen-dependent hypothetical gene (PVX_114580).
Using CRISPR/Cas9, PVX_114580 was integrated into P.falciparum 3D7, transcribed, and translated during asexual stages. Immunofluorescence analysis demonstrated that the protein, proposed to be named P. vivax Spleen-Dependent Protein 1
(PvSDP1), was located at the surface of infected red blood cells and this localization was confirmed in natural infections. Plasma-derived EVs from P. vivax-infected individuals significantly increased cytoadherence of transgenic parasites to hSFs and this binding was inhibited by anti-PvSDP1 antibodies. Single-cell RNAseq of PvEVs-treated hSFs revealed increased expression of adhesion-related genes. These findings demonstrate the importance of parasite spleen-dependent genes and EVs from natural infections in the formation of intrasplenic niches in P. vivax, a major challenge for malaria elimination.

**Single Cell RNA Sequencing from: PvSDP1 and its role in Extracellular Vesicles-Mediated Intrasplenic Infections in Plasmodium vivax**

In this repository the code generated to analyzed the 10X data generated in the paper "PvSDP1 and its role in Extracellular Vesicles-Mediated Intrasplenic Infections in Plasmodium vivax" is upload. 
The analysis where done using: R version 4.2.2 (2022-10-31)
The main library used to analyzed was Seurat Version 4.9.9
A total of 10,338 cells in the EV-stimulated sample and 8,894 cells in the non-stimulated sample where analyzed

**Summary of the pipeline:**
1. Quality control,removing cells with fewer than 300 genes and those with a mitochondrial RNA content exceeding 10%.
2. Doublets were identified using the DoubletFinder library (Version 2.0.3) integrated into the Seurat pipeline, to discern doublets within each sample.
3. Normalization of the dataset using LogNormalized () function
4. Choosing highly variable genes
5. Estimating pca dimmensions
6. Clustering analysis
7. Identify top differentially expressed genes between studied conditions. 
