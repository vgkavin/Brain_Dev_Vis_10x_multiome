library(Seurat)
library(zellkonverter)
library(SeuratDisk)
library(SingleCellExperiment)
library(Signac)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(pheatmap)
library(stringr)
library(EnsDb.Mmusculus.v79)
library(biovizBase)
library(Matrix)


#load the ATACseq sce dataset
vis.atac.sce <- readRDS("DevVIS_multiome_snATAC_processed.rds")


#identify matching cells using barcodes to test multiome validity
#Extract cell barcodes
rna.barcodes <- str_extract(colnames(vis.rna.sce), "^[A-Z]+")
atac.barcodes <- str_extract(colnames(vis.atac.sce), "(?<=#)[A-Z]+")

#Find common barcodes
common.barcodes <- intersect(rna.barcodes, atac.barcodes)
#170k cells share barcodes- Valid for multiome analysis

#subset ~20% of matching cells for analysis to reduce memory usage
selected.barcodes <- sample(common.barcodes, 35000) 

#Subset & Align both objects
vis.rna.sce <- vis.rna.sce[, match(selected.barcodes, rna.barcodes)]
vis.atac.sce <- vis.atac.sce[, match(selected.barcodes, atac.barcodes)]

#snRNAseq data has extensive metadata while ATAC lacks any
#Transfer some vital metadata
target_cols <- c("library_prep", "roi", "sex", "age_label")
colData(vis.atac.sce) <- colData(vis.rna.sce)[, target_cols]


#the ATACseq sce object is still huge as it is a dense matrix and not sparse
#run the atac processing.R script to convert the large dense matrix to a memory efficient sparse matrix


#after processing ATAC seq sce data into sparse matrix
#now single cell experiment objects for both datasets with the same set of 35k cells is available


#QC and analysis for transcriptomics data

#create seurat object of subsetted dataset
vis.rna <- CreateSeuratObject(counts = counts(vis.rna.sce))

#Calculate %mitochondrial genes and add to metadata
vis.rna$MT <- PercentageFeatureSet(vis.rna, pattern = "^MT-") 


#QC
#Filtering sets based on criteria provided by author
vis.rna <- subset(vis.rna, subset = nFeature_RNA >= 300 & MT <= 25)

#normalize data
vis.rna <- NormalizeData(vis.rna)

#find differentally expressed genes to segregate clusters
vis.rna <- FindVariableFeatures(vis.rna, nfeatures = 3000, verbose = TRUE)

#scaling to reduce unwanted noise
vis.genes <- rownames(vis.rna)
vis.rna <- ScaleData(vis.rna, features = vis.genes)

#Linear Dimensionality reduction by identifying Princple components 
vis.DEG <- VariableFeatures(vis.rna)
vis.rna <- RunPCA(vis.rna, features = vis.DEG)
ElbowPlot(vis.rna)

#identify genes with similar expression patterns and make cell clusters at author provided resolution
vis.rna <- FindNeighbors(vis.rna)
vis.rna <- FindClusters(vis.rna, resolution = 0.13)
Idents(vis.rna) <- vis.rna$RNA_snn_res.0.13

#Non-linear dimentionality reduction using UMAP
vis.rna <- RunUMAP(vis.rna, dims = 1:20)

#visualize and save
vis.rna.umap <- DimPlot(vis.rna, reduction = "umap", group.by = "RNA_snn_res.0.13")
ggsave(vis.rna.umap, dpi = 300, height = 6, width = 12, filename = "Dev_Vis_RNA_UMAP.png")

# Identify DEGs for cluster "3" 
vis.rna.DE.c3 <- FindMarkers(
  object = vis.rna,          
  ident.1 = 3,               
  only.pos = FALSE,          
  min.pct = 0.1,             
  logfc.threshold = 1     
)

#save results
write.csv(vis.rna.DE.c3, "DEGs_cluster3.csv")


#QC and analysis for ATAC seq data

#load data if not availabkle in local environment
vis.atac.sce <- readRDS("vis_atac_35k_sparse_final.rds")

#extract peaks data
peak_data <- rowData(vis.atac.sce)


#identify peak ranges for chromatin assay
peak_ranges <- GRanges(
  seqnames = peak_data$seqnames,
  ranges = IRanges(start = peak_data$start, end = peak_data$end), # Assumes 'end' exists
  strand = peak_data$strand
)


#extract ATAC counts data
counts.atac <- counts(vis.atac.sce)
#extract metadata
metadata.atac <- as.data.frame(colData(vis.atac.sce))

#create chromatic assay object with peaks, counts and metadata
vis.atac <- CreateChromatinAssay(
  counts = counts.atac,
  ranges = peak_ranges,
  min.cells = 10,
  min.features = 200
)

#store the chromatin assay data in a seurat object
vis.atac.obj <- CreateSeuratObject(
  counts = vis.atac,
  metadata = metadata.atac,
  assay = "ATAC"
)


#extract annotation data and map it to metadata
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)


seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

vis.atac.obj@assays$ATAC@annotation <- annotation



#QC
#filtering low quality cells based on counts
#could not filter based on fragment data as frag.file or data was not available
vis.atac.obj  <- subset(x = vis.atac.obj, subset = nCount_ATAC > 3000 & nCount_ATAC < 30000 )

#normalization
vis.atac.obj <- RunTFIDF(vis.atac.obj)

#find enriched peaks
vis.atac.obj <- FindTopFeatures(vis.atac.obj, min.cutoff = "q0")

#linear dimensionality reduction
vis.atac.obj <- RunSVD(vis.atac.obj)

#cell clustering
vis.atac.obj <- FindNeighbors(object = vis.atac.obj, reduction = "lsi", dims = 2:30)
vis.atac.obj <- FindClusters(object = vis.atac.obj, algorithm = 3, resolution = 0.13)

#non linear dimensionality reduction
vis.atac.obj <- RunUMAP(vis.atac.obj, reduction = "lsi", dims = 2:30)

#visualize and save
vis.atac.umap <- DimPlot(vis.atac.obj, group.by = "ATAC_snn_res.0.13")
ggsave(vis.atac.umap, dpi = 300, height = 6, width = 12, filename = "Dev_Vis_ATAC_UMAP.png")


vis.atac.chrom.c3 <-FindMarkers(
  object = vis.atac.obj,
  ident.1 = 3,
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  min.pct = 0.1,
  logfc.threshold = 1,
  only.pos = TRUE
  )

write.csv(vis.atac.chrom.c3, "enriched-open_chromatin_cluster3.csv")

saveRDS(vis.atac.obj, "vis.seurat.rds")
vis.atac.obj <- readRDS("vis.seurat.rds")
