#library(BiocManager)
#BiocManager::install("ComplexHeatmap")
remotes::install_github("bioFAM/MOFA2", ref = "e605bf2", upgrade = "never")
remotes::install_github("saezlab/MOFAcellulaR", ref = "8cb0785", upgrade = "never")

library(MOFA2)
library(MOFAcellulaR)
