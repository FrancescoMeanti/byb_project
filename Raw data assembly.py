import pandas as pd
import scanpy as sc


adata = sc.read_mtx(r"data\raw data\RNA_counts.mtx")
adata_bc=pd.read_csv(r"data\raw data\cells.csv",header=0)
adata_features=pd.read_csv(r"data\raw data\genes.csv",header=0)
adata= adata.T
adata.obs['cell_id']= adata_bc['cell'].tolist()
adata.var['gene_name']= adata_features['gene'].tolist()
adata.var_names = adata_features['gene'].tolist()
sample_obs = pd.read_csv(r"data\raw data\Combined_Metadata.csv",header=0)
adata.obs = sample_obs
adata.raw = adata.copy()
 
 
#adata.var_names_make_unique()
immune_data = adata

immune_data.write_h5ad(r"data\H5AD Objects\Immune_Raw_New.h5ad")
