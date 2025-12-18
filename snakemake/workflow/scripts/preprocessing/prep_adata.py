import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

input = snakemake.input[0]
study = snakemake.wildcards['study']
cmap_cells = snakemake.params[0]
cmap_conditions = snakemake.params[1]
output = snakemake.output[0]
output_pdf =snakemake.output[1]

adata =  sc.read_h5ad(input)
adata.obs['annotation_MOFA'] = None
adata.obs['cond_test']='control'


# each study is in a slightluy separate format, the following lines unify formats
if study == 'schiller':
    adata.obs.loc[adata.obs.louvain_krt.str.contains('fibroblast', case = False),'annotation_MOFA'] = 'fibroblast'
    adata.obs.loc[adata.obs.louvain_krt.str.contains('Krt|AT1|AT2', case = False),'annotation_MOFA'] = 'epithelial'
    adata.obs.loc[adata.obs.louvain_krt.str.contains('mono|macro|Am|dendritic', case = False),'annotation_MOFA'] = 'monomacro'
    adata.obs.loc[adata.obs.louvain_krt.str.contains('endoth', case = False),'annotation_MOFA'] = 'endothelial'
    adata.obs.loc[adata.obs.grouping!='PBS','cond_test'] = 'fibrosis'
    adata.obs['sample'] = adata.obs['identifier'] 
    adata_new = adata[adata.obs['annotation_MOFA'].notna()].copy()
    adata_new.obs['study'] = 'schiller'
    

    
elif study == 'peyser':
    adata.obs.loc[adata.obs.annotation.str.contains('fibroblast', case = False),'annotation_MOFA'] = 'fibroblast'
    adata.obs.loc[adata.obs.annotation.str.contains('AT2|AT1|Krt8+', case = False),'annotation_MOFA'] = 'epithelial'
    adata.obs.loc[adata.obs.annotation_macro.str.contains('Md|Macro|Mono', case = False),'annotation_MOFA'] = 'monomacro'
    adata.obs.loc[adata.obs.annotation.str.contains('endoth', case = False),'annotation_MOFA'] = 'endothelial'
    adata.obs.loc[adata.obs.cond!='1','cond_test'] = 'fibrosis'
    adata.obs['sample'] = adata.obs['batch'] 
    adata_new = adata[adata.obs['annotation_MOFA'].notna()].copy()
    adata_new.obs['study'] = 'peyser'
    
    
    
elif study == 'xie':
    adata.obs.loc[adata.obs.louvain_final.str.contains('fibro', case = False),'annotation_MOFA'] = 'fibroblast'
    adata.obs.loc[adata.obs.louvain_macro2.str.contains('macro|Am|mono|dend', case = False),'annotation_MOFA'] = 'monomacro'
    adata.obs.loc[adata.obs.louvain_final.str.contains('endoth', case = False),'annotation_MOFA'] = 'endothelial'
    adata.obs['cond_test']='control'
    adata.obs.loc[adata.obs.cond!='control','cond_test'] = 'fibrosis'
    adata.obs['study'] = 'xie'
    adata.obs['sample'] = 'sample_' + adata.obs['sample'].astype(str)    
    adata_new = adata[adata.obs['annotation_MOFA'].notna()].copy()
    

    
elif study == 'tsukui':
    adata.obs.loc[adata.obs.louvain_myofibro.str.contains('fibro', case = False),'annotation_MOFA'] = 'fibroblast'
    adata.obs.loc[adata.obs.louvain_myofibro.str.contains('Krt8|AT1|AT2', case = False),'annotation_MOFA'] = 'epithelial'
    adata.obs.loc[adata.obs.louvain.str.contains('mono|macro', case = False),'annotation_MOFA'] = 'monomacro'
    adata.obs.loc[adata.obs.louvain.str.contains('endoth', case = False),'annotation_MOFA'] = 'endothelial'
    adata.obs.loc[adata.obs.cond=='0', 'cond_test']='fibrosis'
    adata.obs['sample'] = [i[-7:-5] for i in adata.obs['batch']]
    adata_new = adata[adata.obs['annotation_MOFA'].notna()].copy()
    adata_new.obs['study'] = 'tsukui'
    
    
    
elif study == 'misharin':
    adata.obs.loc[adata.obs.louvain_fibro.str.contains('fibro', case = False),'annotation_MOFA'] = 'fibroblast'
    adata.obs.loc[adata.obs.louvain_fibro.str.contains('Krt|AT', case = False),'annotation_MOFA'] = 'epithelial'
    adata.obs.loc[adata.obs.louvain_fibro.str.contains('Mono|macro|Am|dend|pDC', case = False),'annotation_MOFA'] = 'monomacro'
    adata.obs.loc[adata.obs.louvain_fibro.str.contains('endoth', case = False),'annotation_MOFA'] = 'endothelial'
    adata.obs.loc[adata.obs.batch=='1', 'cond_test']='fibrosis'
    adata.obs['sample'] = adata.obs['cond']
    adata_new = adata[adata.obs['annotation_MOFA'].notna()].copy()
    adata_new.obs['study'] = 'misharin'
     
    
else:
    raise ValueError('submitted file does not have a specification how to be processed')


    
with PdfPages(output_pdf) as output_pdf:
    for i in [adata, adata_new]:
        fig, axs = plt.subplots(1,3, tight_layout = True, figsize = (15,4))
        sc.pl.umap(i, color = ['annotation_MOFA'], ax = axs[0], show = False, palette = cmap_cells)
        sc.pl.umap(i, color = ['cond_test'], ax = axs[1], show = False, palette = cmap_conditions)
        sc.pl.umap(i, color = ['sample'], ax = axs[2], show = False)
        output_pdf.savefig(fig)
        
        
adata_new.write(output)

