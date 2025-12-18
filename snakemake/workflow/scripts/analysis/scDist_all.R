library(scDist)
library(anndata)
library(magrittr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
set.seed(1126490984)

# Enable parallelization
plan(sequential)
options(future.globals.maxSize= 200000000000)

# snakemake input and output
h5ad_file = snakemake@input$adata_files
results_path = snakemake@output$results
feature_importance_path = snakemake@output$feature_importance

file <- h5ad_file

# Step 1: Read the h5ad file
ad <- read_h5ad(file)

# Step 2: Filter the rows where 'annotation_MOFA' is not NA
ad <- ad[!is.na(ad$obs$annotation_MOFA), ]

# filter cell types that do not contain cells from disease & reference
fib_ctypes <- unique(ad$obs[(ad$obs$cond_test == 'fibrosis'), ]$annotation_MOFA)
ref_ctypes <- unique(ad$obs[(ad$obs$cond_test == 'control'), ]$annotation_MOFA)

ctypes <- as.character(unique(fib_ctypes[fib_ctypes %in% ref_ctypes]))
ad <- ad[(ad$obs$annotation_MOFA  %in% ctypes), ]



# Step 3: Extract metadata and counts
metadata <- rbind(ad$obs)[c('sample', 'cond_test', 'annotation_MOFA', 'study')]
Y <- ad$layers["counts"]



print(head(metadata))


# Step 4: Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = t(Y), meta.data = metadata)

# Clean up memory
rm(ad, metadata, Y)
gc()

# Step 5: Apply SCTransform normalization
seurat_obj <- SCTransform(seurat_obj, ncells=5000, conserve.memory=TRUE)


# run scDist
out <- scDist(seurat_obj@assays$SCT@scale.data,seurat_obj@meta.data,fixed.effects = "cond_test",
              random.effects=c("sample"),
              clusters="annotation_MOFA")


# access scDist results
scd.object <- out
results <- scd.object$results



# summary
result_tibble <- as_tibble(lapply(ctypes, function(name) scd.object$vals[[name]]$beta.hat), .name_repair  = 'unique')
# set column names
names(result_tibble) <- ctypes
result_tibble <- result_tibble %>%
  mutate("gene" = scd.object$gene.names)

# save results
write.table(result_tibble , file = feature_importance_path)
write.table(results , file = results_path)
