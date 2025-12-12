#initial script to load the .h5ad anndata object onto a python  virtual env,
#manipulate dataset to be compatible with seurat's format
#the data has to be transposed or any functions like Convert from SeuratDisk or Zellkonvert will throw an error or a fully zero counts matrix

#load required libraries
library(reticulate)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)



#Create a NEW, clean py virtual environment

virtualenv_create("clean-analysis-env")

#Install packages (Force NumPy 1.23.5 for stability)
virtualenv_install("clean-analysis-env", 
                   packages = c("anndata", "numpy==1.23.5", "scipy", "pandas"))

#Activate the new environment
use_virtualenv("clean-analysis-env", required = TRUE)


#verifying the env config
py_config()


#load the anndata object in python and convert counts to doubles, extract counts data and metadata for transposing
# We write a Python script as a string and run it to ensure it is not run on R env
py_run_string("
import anndata
import scipy.sparse
import numpy as np

adata = anndata.read_h5ad('DevVIS_multiome_snRNA_processed.h5ad')

X = adata.X.astype(np.float64)
coo = X.tocoo()
row = coo.row
col = coo.col
data = coo.data

# Extract names and metadata
obs_names = adata.obs_names.to_list()
var_names = adata.var_names.to_list()
meta_data = adata.obs
shape = adata.shape
")


#QC to validate py script execution
cat("Row length:", length(py$row), "\n")
cat("Col length:", length(py$col), "\n")
cat("Data length:", length(py$data), "\n")


#recreate the indices to be 100% sure they are correct.
row_indices <- as.numeric(py$row) + 1
col_indices <- as.numeric(py$col) + 1
data_values <- as.numeric(py$data)

# Get dimensions 
n_rows_py <- as.integer(py$shape[[1]]) # Python rows (Cells)
n_cols_py <- as.integer(py$shape[[2]]) # Python cols (Genes)

#make a counts object with extracted counts and dims data
counts <- sparseMatrix(
  i = as.numeric(py$row) + 1,
  j = as.numeric(py$col) + 1,
  x = as.numeric(py$data),
  dims = c(n_rows_py, n_cols_py)
)


# transpose
counts <- Matrix::t(counts)


#assign names for each part of seurat object
cell_names <- py$obs_names
gene_names <- py$var_names

colnames(counts) <- cell_names
rownames(counts) <- gene_names

meta_data <- py$meta_data

#create seurat object with wrangled data
vis.rna <- CreateSeuratObject(counts = counts, meta.data = meta_data)

#validate the seurat obj
View(vis.rna@meta.data)
View(vis.rna$cellNames)

#convert to single cell experiment to match with ATAC seq data and subset
vis.rna.sce <- as.SingleCellExperiment(vis.rna)

#continue analysis in analysis.R file
