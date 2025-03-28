import scanpy as sc
import pandas as pd

immune_data = sc.read_h5ad(r"data\H5AD Objects\Immune_Raw_New_Harmonized.h5ad")
immune_data = immune_data[immune_data.obs["cell_type"] != 'Erythrocyte']

# Filter to only Early and Late groups
immune_data = immune_data[immune_data.obs["treatment_group"].isin(["WoundEarly", "WoundLate"])].copy()
ensembl = pd.read_csv("Ensembl_New.txt")

# Save the original ENSEMBL IDs
original_ids = immune_data.var_names.copy()
# Create dictionaries containing salmon, zebrafish and human gene names
salmon_to_gene = dict(zip(ensembl["Gene stable ID"], ensembl["Gene name"]))
zebra_to_gene = dict(zip(ensembl["Gene stable ID"], ensembl["Zebrafish gene name"]))
human_to_gene = dict(zip(ensembl["Gene stable ID"], ensembl["Human gene name"]))
var_names_series = pd.Series(immune_data.var_names, index=immune_data.var_names)

# Map the dictionaries to the Ensembl IDs stored in .var
immune_data.var["salmon_name"] = original_ids.map(salmon_to_gene)
immune_data.var["zebra_name"] = original_ids.map(zebra_to_gene)
immune_data.var["human_name"] = original_ids.map(human_to_gene)
immune_data.var["final_name"] = immune_data.var["salmon_name"].fillna(immune_data.var["zebra_name"]).fillna(immune_data.var["human_name"]).fillna(var_names_series)
# Create a final mapping dictionary
gene_id_to_name = dict(zip(original_ids, immune_data.var["final_name"]))


# Perform differential expression analysis
sc.tl.rank_genes_groups(immune_data, groupby="treatment_group", method="wilcoxon")

clusters = immune_data.obs["treatment_group"].unique()

output_file = "Ranked genes by TG.xlsx"

 # Define the desired order of columns
desired_order = ["Ensembl ID", "Gene Name", "logfoldchanges", "scores", "pvals", "pvals_adj"]

with pd.ExcelWriter(output_file, engine="xlsxwriter") as writer:
    # Extract rank_genes_groups data
    for cluster in clusters:
        ranked_genes = sc.get.rank_genes_groups_df(immune_data, group = cluster, pval_cutoff= 0.0005)
        ranked_genes["Gene Name"] = ranked_genes["names"].map(gene_id_to_name)
        ranked_genes.rename(columns={'names': 'Ensembl ID'}, inplace=True)
        # Reorder the DataFrame columns
        ranked_genes = ranked_genes[desired_order]

        ranked_genes.to_excel(writer, sheet_name=f"{cluster}", index=False)