# 🧬 Shared and organ-specific gene expression signatures of fibrotic diseases
## 🔍 Overview
<img width="300" height="60" alt="organs" src="https://github.com/user-attachments/assets/22860868-1c0c-47cf-987e-b5edc295d0c4" />

This repository contains the code corresponding to our manuscript:
> *citation for preprint*

We conducted a large-scale meta-analysis of single-cell transcriptomic data from human healthy and fibrotic tissues to identify both shared and organ-specific transcriptomic profiles. Using datasets from the heart, kidney, lung, and liver, we constructed a single-cell fibrosis atlas of over five million cells from 20 studies, covering more than 25 etiologies across four organs.  
Through systematic comparison of these datasets, we identified organ-specific as well as cross-organ fibrotic gene expression profiles in major cell types and disease fibroblast subpopulations, characterized by the excessive production of extracellular matrix, revealing a shared fibrotic response across tissues. 

## 💻 Code overview
The code in structured into a snakamake pipeline. Configurations are found in the `/profile` directory. Snakemake rules, envrionments and analysis code is found in the `/workflow` directory. 

The code was modulatized into the following modules:

- preprocessing  
    - processing & harmonization of raw datasets 
- analysis
    - core analysis of scRNA/snRNA datasets
- integration  
    - GPU-dependent steps (for cluster configurations, see `/profile/slurm2` - contains cell type annotation transfer and mesenchymal cell integration  
- myofib  
     - processing & analysis of disease fibroblast subset of the data
- spatial
     - preprocessing & analysis of spatial datasets
- plotting
     - final processing and visualization of data
- dataformat
     - formatting of data for publication


## 📂 Where to get the data
Non-processed data is found in the publications of indiviual datasets. 
Processed pseudobulks and analysis results can be accessed in zenodo: 
> *zenodo link*

## 🌐 User-interface to explore the data
Have a look at our interactive website where you can look up how your gene of interest behaves in fibrotic disease tissues:
> *website link*
