install.packages("anndata", repos = "http://cran.us.r-project.org")

library(BiocManager)
BiocManager::install('zellkonverter')
BiocManager::install('miloR')
BiocManager::install('SingleCellExperiment')
BiocManager::install('scater')
BiocManager::install('scran')
BiocManager::install("S4Vectors")



library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(zellkonverter)

