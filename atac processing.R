library(Matrix)
library(SingleCellExperiment)
library(parallel)

#continued from analysis.R after matching and subsetting 35k cells in both datasets
#the ATACseq sce object is still huge as it is a dense matrix and not sparse
#convert it to sparse to save momory and speed up analysis
#spliting 35k cells into chunks of 2k cells each and converting them into sparse matrices simultaneously on 8 cores
cores <- 8 
chunk_size <- 2000
total_cells <- ncol(vis.atac.sce)

# Define the chunks
chunk_indices <- split(seq_len(total_cells), ceiling(seq_len(total_cells)/chunk_size))
num_chunks <- length(chunk_indices)

#Define the conversion function
# This function will run on 8 cores simultaneously
convert_chunk <- function(indices, mat_obj) {
  # Extract the slice
  slice <- assay(mat_obj, "counts")[, indices]
  # Convert to sparse
  as(slice, "dgCMatrix")
}

# run the function for parallel processing of 8 chunks at a time
sparse_chunks <- mclapply(chunk_indices, 
                          convert_chunk, 
                          mat_obj = vis.atac.sce, 
                          mc.cores = cores)

#Combine and save
final_sparse <- do.call(cbind, sparse_chunks)

# Replace the heavy matrix with sparse data
assay(vis.atac.sce, "counts") <- final_sparse

#clear memorny usage by environment
rm(final_sparse, sparse_chunks)
gc()

#save RDS for further processing
saveRDS(vis.atac.sce, "vis_atac_35k_sparse_final.rds")

#continue downstream analysis in analysis.R