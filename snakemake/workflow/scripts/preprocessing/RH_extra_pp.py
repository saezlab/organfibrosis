# script to preprocess Reichart data


import scanpy as sc
import numpy as np
import pandas as pd


input_file = snakemake.input[0]


output = snakemake.output[0]


sc_dat = sc.read_h5ad(filename = input_file)
sc_dat.obs["barcode"] = sc_dat.obs.index.values
sc_dat.X = sc_dat.raw.X
# keep only assigned cells
sc_dat = sc_dat[sc_dat.obs[["Assigned"]].values == True,]
disease_df = {'disease': ['dilated cardiomyopathy', 
                          'normal',
                          'arrhythmogenic right ventricular cardiomyopathy',
                          'non-compaction cardiomyopathy'],
              'disease_code': ["DCM","NF","ARVC","NCC"],
              'heart_failure': ["HF", "NF", "HF", "HF"]}

disease_df = pd.DataFrame(disease_df)

new_codes = sc_dat.obs.merge(disease_df, on='disease', 
                           how='left')[["barcode","disease_code","heart_failure"]]
new_codes.set_index("barcode", inplace = True)
new_codes = new_codes.loc[sc_dat.obs.index.values, :]

sc_dat.obs["disease_code"] = new_codes["disease_code"].values
sc_dat.obs["heart_failure"] = new_codes["heart_failure"].values

# Now to make the dictionary of matched cell types
available_cells = sc_dat.obs[["cell_type"]].drop_duplicates()
available_cells["cell_type_uni"] = "none"

# assign cell types
available_cells.loc[available_cells["cell_type"].str.contains('muscle'), 'cell_type_uni'] = "CM"
available_cells.loc[available_cells["cell_type"].str.contains('endothelial'), 'cell_type_uni'] = "Endo"
available_cells.loc[available_cells["cell_type"].str.contains('fibroblast'), 'cell_type_uni'] = "Fib"
available_cells.loc[available_cells["cell_type"].str.contains('lymphocyte'), 'cell_type_uni'] = "Lymphoid"
available_cells.loc[available_cells["cell_type"].str.contains('mast'), 'cell_type_uni'] = "Myeloid"
available_cells.loc[available_cells["cell_type"].str.contains('myeloid'), 'cell_type_uni'] = "Myeloid"


new_cts = sc_dat.obs.merge(available_cells, on='cell_type', 
                           how='left')[["barcode","cell_type_uni"]]
new_cts.set_index("barcode", inplace = True)
new_cts = new_cts.loc[sc_dat.obs.index.values, :]

sc_dat.obs["cell_type_uni"] = new_cts["cell_type_uni"].values
sc_dat.obs.loc[sc_dat.obs["cell_states"].str.contains('PC'), 'cell_type_uni'] = "PC"
sc_dat.obs.loc[sc_dat.obs["cell_states"].str.contains('SMC'), 'cell_type_uni'] = "vSMCs"
sc_dat = sc_dat[sc_dat.obs[["cell_type_uni"]].values != "none",]

sc_dat.obs['tissue'] = sc_dat.obs['tissue'].cat.rename_categories({
    'heart left ventricle':'LV', 
    'heart right ventricle':'RV', 
    'apex of heart':'A', 
    'interventricular septum':'S'})


# Filter obs to contain things that are relevant
sc_dat.obs = sc_dat.obs[['Sample',
                         'Primary.Genetic.Diagnosis',
                         'cell_type_uni',
                         'assay',
                         'suspension_type',
                         'disease',
                         'sex',
                         'disease_code',
                         'heart_failure',
                         'tissue',
                         'donor_id',
                        'cell_states']]

# rename columns to fit other studies
sc_dat.obs = sc_dat.obs.rename(columns={"Sample": "sample", 
                                        "cell_type_uni":"cell_type", 
                                        'tissue':'region', 
                                        'assay':'tech',
                                        'cell_states':'cell_type2'})


sc_dat.write(output)