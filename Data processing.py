import scanpy as sc
import scanpy.external as sce

sc.settings.set_figure_params(dpi=50, facecolor="white")


immune_data = sc.read_h5ad("Immune_Raw_New.h5ad")


sc.pp.filter_cells(immune_data, min_genes=100)
sc.pp.filter_genes(immune_data, min_cells=3)


# Save raw count data to a new layer called "raw_counts"
immune_data.layers["raw_counts"] = immune_data.X.copy()


# Normalizing to median total counts
sc.pp.normalize_total(immune_data)
# Logarithmize the data
sc.pp.log1p(immune_data)

sc.tl.pca(immune_data)

# Run Harmony
sce.pp.harmony_integrate(
    immune_data,
    "ident",
    basis = "X_pca",
    adjusted_basis= "X_pca_harmony" # Field in which adjusted PCA will be stored
)

# Compute neighbors using corrected PCA
sc.pp.neighbors(immune_data, use_rep="X_pca_harmony")

# Compute UMAP
sc.tl.umap(immune_data)

# Aplly Leiden clustering to cells
sc.tl.leiden(immune_data, key_added="leiden", resolution=0.07, flavor="igraph")


# Map cluster names to leiden clusters
cluster_rename_dict = {
    '0': 'MITFa+ Phagocyte',
    '1': 'T Lymphocyte',
    '2': 'Mononuclear Phagocyte',
    '3': 'Dendritic Cell-like',
    '4': 'B Lymphocyte',
    '5': 'Erythrocyte',
    '6': 'Granulocyte',
}

immune_data.obs['cell_type'] = immune_data.obs['leiden'].map(cluster_rename_dict).astype('category')

immune_data.write_h5ad("Immune_Raw_New_Harmonized.h5ad")
