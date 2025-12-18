install.packages("anndata", repos = "http://cran.us.r-project.org")

library(BiocManager)
BiocManager::install('glmGamPoi')

remotes::install_github("bioFAM/MOFA2", ref = "e605bf2", upgrade = "never")
devtools::install_github("phillipnicol/scDist", ref = "4b54dd7", upgrade = "never")

library(anndata)
library(scDist)
